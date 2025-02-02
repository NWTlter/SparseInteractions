# This script will fit the model with a monotonic lambda and constant alphas
#    (with respect to the environment) to the corresponding data, using a series
#    of dataset sizes (10, 20, 50, 100, 200), and then examine the accuracy of
#    parameter estimates and posterior predictive checks on 300 out of sample
#    data points.

rm(list = ls())
library(here)
library(rstan)
library(HDInterval)
library(RColorBrewer)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the current sample size and associated prefix for all graph and result
#    file names
N <- 200
max_N <- 200
FilePrefix <- paste("N", N, "_", sep = "")

# Now assign the focal species and the file paths for the stan models
# These paths are within the "Box" folder (this file's location), and may need to be updated
# to the user's file structure
Focal <- 8
PrelimStanPath <- here("BH_simulations/Box/StanCode/Prelim_monoLambda_constAlpha.stan")
FinalStanPath <- here("BH_simulations/Box/StanCode/Final_monoLambda_constAlpha.stan")

# Load in the appropriate data
FullSim <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/simulation_const.csv")) %>%
        select(-X)
TrueVals <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/parameters_const.csv"))
TrueAlphas <- TrueVals$alpha.1

# assign some universal values to be used across model fits and graphs
S <- 15
Intra <- rep(0, S)
Intra[Focal] <- 1
tau0 <- 1
slab_df <- 4
slab_scale <- sqrt(2)

# Set initial values to avoid initial problems with the random number generator
ChainInitials <- list(lambdas = c(TrueVals$lambda.mean[Focal], TrueVals$lambda.env[Focal]), 
                      alpha_generic_tilde = mean(TrueAlphas)/0.75 + 1,  
                      alpha_hat_ij_tilde = rep(0, S), c2_tilde = 1.25,
                      local_shrinkage_ij = rep(5, S), tau_tilde = 15,
                      alpha_intra_tilde = TrueAlphas[Focal])
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# save the values for the posterior predictive checks
ppc_data <- subset(FullSim, (species == Focal) & (run > max_N) & (time == 0) & (thinned == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
     SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) &
                                        (time == 0) & (thinned == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

# Create the data vectors to be passed to rstan for subsequent model fits
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij")

# Set the local values to pass to rstan
FullData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0) & (thinned == 0))
ThinData <- subset(FullSim, (species == Focal) & (run <= N) & (time == 0) & (thinned == 1))
Nt <- c(FullData$pop, ThinData$pop)
env <- c(FullData$run.env, ThinData$run.env)
SpMatrix <- matrix(data = NA, nrow = 2*N, ncol = S)
for(s in 1:S){
        SpMatrix[1:N,s] <- subset(FullSim, (species == s) & (run <= N) & (time == 0) & (thinned == 0))$pop
        SpMatrix[(N+1):(2*N),s] <- subset(FullSim, (species == s) & (run <= N) & (time == 0) & (thinned == 1))$pop
}
Ntp1 <- c(subset(FullSim, (species == Focal) & (run <= N) & (time == 1) & (thinned == 0))$pop,
          subset(FullSim, (species == Focal) & (run <= N) & (time == 1) & (thinned == 1))$pop)

# Now run the perliminary fit of the model to assess parameter shrinkage
N <- length(Nt)
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, init = InitVals, control = list(max_treedepth = 15, adapt_delta = 0.995))
PrelimPosteriors <- extract(PrelimFit)
FitFileName <- paste(here("BH_simulations/Box/StanFits/monoLambda_constAlpha/"), FilePrefix, "PrelimFit.rdata", sep = "")
save(PrelimFit, PrelimPosteriors, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
traceplot(PrelimFit, pars = "lambdas")
traceplot(PrelimFit, pars = "alpha_generic")
traceplot(PrelimFit, pars = "alpha_intra")
traceplot(PrelimFit, pars = "alpha_hat_ij")
acf(PrelimPosteriors$alpha_generic)
acf(PrelimPosteriors$alpha_intra)
PlotSamples <- sample(1:S, size = 4, replace = FALSE)
for(i in 1:4){
        acf(PrelimPosteriors$alpha_hat_ij[,PlotSamples[i]])
}

### Final model -----------
# Determine the parameters that should be included and run the final model
plot(PrelimFit, pars = "alpha_hat_ij")

Inclusion_ij <- rep(0, S)
IntLevel <- 0.5
for(s in 1:S){
     Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
     if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
          Inclusion_ij[s] <- 1
     }
     if(s == Focal){
             Inclusion_ij[s] <- 0
     }
}
Inclusion_ij

# Reset initial conditions to allow faster fitting
ChainInitials <- list(lambdas = colMeans(PrelimPosteriors$lambdas), 
                      alpha_generic_tilde = mean(PrelimPosteriors$alpha_generic_tilde), 
                      alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), 
                      local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
                      c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
                      alpha_intra_tilde = mean(PrelimPosteriors$alpha_intra_tilde))
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals, control = list(max_treedepth = 15))
FinalPosteriors <- rstan::extract(FinalFit)
FitFileName <- paste(here("BH_simulations/Box/StanFits/monoLambda_constAlpha/"), FilePrefix, "FinalFit.rdata", sep = "")
save(FinalFit, FinalPosteriors, Inclusion_ij, file = FitFileName)

# Examine diagnostics and determine if parameters of model run should be updated
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = "lambdas")
traceplot(FinalFit, pars = "alpha_generic")
traceplot(FinalFit, pars = "alpha_intra")
which(Inclusion_ij == 1)
traceplot(FinalFit, pars = "alpha_hat_ij")

# Double check the autocorrelation
acf(FinalPosteriors$lambdas[,1])
acf(FinalPosteriors$lambdas[,2])
acf(PrelimPosteriors$alpha_generic)
acf(PrelimPosteriors$alpha_intra)
for(s in 1:S){
        if(Inclusion_ij[s] == 1){
                acf(FinalPosteriors$alpha_hat_ij[,s])
        }
}


########### Posterior Predictive Check
PostLength <- length(FinalPosteriors$alpha_generic)
# calculate the posterior distributions of the interaction coefficients
alpha_ij <- matrix(NA, nrow = PostLength, ncol = S)
for(s in 1:S){
     if(s == Focal){
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_intra)
     }else{
          alpha_ij[,s] <- exp(FinalPosteriors$alpha_generic + 
                                      Inclusion_ij[s] * FinalPosteriors$alpha_hat_ij[,s])
     }
     
}
alpha_intra <- exp(FinalPosteriors$alpha_intra)
# calculate the posterior distributions of lambda_ei
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
     lambda_ei[,i] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[i])
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_ij[i,] * SpMatrix_ppc[j,])
          Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
          Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
     }
}

# Calculate final fit results for the ppc
PredVals <- matrix(data = NA, nrow = 3, ncol = N_ppc)
for(i in 1:N_ppc){
        PredVals[1,i] <- mean(Growth_pred[,i], na.rm = TRUE)
        PredVals[2:3,i] <- HDInterval::hdi(Growth_pred[,i])
}

# Calculate the accuracy of parameter estimates from the preliminary fits
LambdaEsts <- matrix(data = NA, nrow = 3, ncol = 2)
LambdaEsts[1,1] <- mean(FinalPosteriors$lambdas[,1] - TrueVals$lambda.mean[Focal])
LambdaEsts[2:3,1] <- hdi(FinalPosteriors$lambdas[,1] - TrueVals$lambda.mean[Focal])
LambdaEsts[1,2] <- mean(FinalPosteriors$lambdas[,2] - TrueVals$lambda.env[Focal])
LambdaEsts[2:3,2] <- hdi(FinalPosteriors$lambdas[,2] - TrueVals$lambda.env[Focal])

# Now calculate the alpha estimates
AlphaEsts <- matrix(data = NA, nrow = 3, ncol = S)
for(s in 1:S){
        if(s == Focal){
                # First the intercept
                Intercept <- FinalPosteriors$alpha_intra
                AlphaEsts[1,s] <- mean(Intercept - TrueAlphas[s])
                AlphaEsts[2:3,s] <- hdi(Intercept - TrueAlphas[s])
        }else{
                # First the intercept
                Intercept <- FinalPosteriors$alpha_generic + FinalPosteriors$alpha_hat_ij[,s] * Inclusion_ij[s]
                AlphaEsts[1,s] <- mean(Intercept - TrueAlphas[s])
                AlphaEsts[2:3,s] <- hdi(Intercept - TrueAlphas[s])
        }
}

# Finally, calculate the "true" generic alpha that the model is attempting to estimate
GenericIntercepts <- 1 - Inclusion_ij
GenericIntercepts[Focal] <- 0
TrueGenericIntercept <- log(sum(colSums(SpMatrix) * GenericIntercepts * exp(TrueAlphas)) / sum(colSums(SpMatrix) * GenericIntercepts))

# Finally, save all the necessary results for the figures
FileName <- paste(here("BH_simulations/Box/StanFits/monoLambda_constAlpha/"), FilePrefix, "GraphStuff.rdata", sep = "")
save(PredVals, Growth_ppc, LambdaEsts, AlphaEsts, Inclusion_ij,
     TrueGenericIntercept, 
     file = FileName)