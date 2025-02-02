# This script will fit the model with a monotonic lambda and constant alphas
#    (with respect to the environment) to the corresponding data, using a series
#    of dataset sizes (10, 50, 200), calculate the posterior predictive estimates
#    on 300 out of sample data points, calculate the deviance in parameter 
#    estimates from the true values, and save it all for later plotting.

rm(list = ls())
library(here)
library(HDInterval)
library(RColorBrewer)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set the current sample size (10, 50, or 200) and associated prefix for all graph and result
#    file names
N <- 200
max_N <- 200
FilePrefix <- paste("N", N, "_", sep = "")

# Now assign the file paths for the stan models
PrelimStanPath <- here("BH_simulations/Main/StanCode/Prelim_monoLambda_envAlpha.stan")
FinalStanPath <- here("BH_simulations/Main/StanCode/Final_monoLambda_envAlpha.stan")

# Load in the appropriate data
load(here("BH_simulations/test_multiple_simulations.RData"))
ExampleSim <- simulations[[8]]
TrueVals <- ExampleSim[[1]]
FullSim <- ExampleSim[[2]]
Focal <- which(TrueVals$focal == 1)
TrueAlphaMeans <- TrueVals$alpha.8 #This simulation has species 8 as the focal
TrueAlphaSlopes <- TrueVals$alpha.env 

# assign some universal values to be used across model fits and graphs
S <- 15
Intra <- rep(0, S)
Intra[Focal] <- 1
tau0 <- 1
slab_df <- 4 
slab_scale <- sqrt(2) 

# save the values for the posterior predictive checks
ppc_data <- subset(FullSim, (species == Focal) & (run > max_N) & (time == 0) & (thinned == 0))
ppc_points <- which(ppc_data$pop > 0)
ppc_runs <- ppc_data$run[ppc_points]
N_ppc <- length(ppc_points)
Nt_ppc <- ppc_data$pop[ppc_points]
env_ppc <- ppc_data$run.env[ppc_points]
SpMatrix_ppc <- matrix(data = NA, nrow = N_ppc, ncol = S)
for(s in 1:S){
     SpMatrix_ppc[,s] <- subset(FullSim, (species == s) & (run %in% ppc_runs) & (time == 0) & (thinned == 0))$pop
}
Ntp1_ppc <- subset(FullSim, (species == Focal) & (run %in% ppc_runs) & (time == 1) & (thinned == 0))$pop
Growth_ppc <- log((Ntp1_ppc + 1)/Nt_ppc)

# Create the data vectors to be passed to rstan for subsequent model fits
PrelimDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
FinalDataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij", "Inclusion_eij")

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

# Now run the preliminary fit of the model to assess parameter shrinkage
N <- 2*N
PrelimFit <- stan(file = PrelimStanPath, data = PrelimDataVec, iter = 3000,
                  chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 15))
PrelimPosteriors <- extract(PrelimFit)

####
save(PrelimFit, PrelimPosteriors, file = "~/Desktop/CurFit.rdata")
load("~/Desktop/CurFit.rdata")
####

# Examine diagnostic plots and determine if the model fit is adequate to move
#       forward with the final fit
# First examine the distribution of Rhat and effective sample size values
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
# Visually examine the traceplots for key parameters
traceplot(PrelimFit, pars = "lambdas")
traceplot(PrelimFit, pars = "alpha_generic")
traceplot(PrelimFit, pars = "alpha_intra")
traceplot(PrelimFit, pars = "alpha_hat_ij")
traceplot(PrelimFit, pars = "alpha_hat_eij")
# Check for parameter correlations and divergent transitions
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
# Check for autocorrelation in key parameters
acf(PrelimPosteriors$alpha_generic[,1])
acf(PrelimPosteriors$alpha_generic[,2])
acf(PrelimPosteriors$alpha_intra[,1])
acf(PrelimPosteriors$alpha_intra[,2])
PlotSamples <- sample(1:S, size = 4, replace = FALSE)
quartz()
par(mfrow = c(2,2))
for(i in 1:4){
        acf(PrelimPosteriors$alpha_hat_ij[,PlotSamples[i]])
}
for(i in 1:4){
        acf(PrelimPosteriors$alpha_hat_eij[,PlotSamples[i]])
}

#### If the diagnostic plots don't reveal any problems wiht the model fit, now
#       move on to determining which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
Inclusion_ij <- rep(0, S)
Inclusion_eij <- rep(0, S)
IntLevel <- 0.5
for(s in 1:S){
     Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
     Ints_eij <- hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = IntLevel)
     if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
          Inclusion_ij[s] <- 1
     }
     if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
             Inclusion_eij[s] <- 1
     }
     if(s == Focal){
             Inclusion_ij[s] <- 0
             Inclusion_eij[s] <- 0
     }
}
Inclusion_ij
Inclusion_eij

# Reset initial conditions with values from the preliminary fit
ChainInitials <- list(lambdas = colMeans(PrelimPosteriors$lambdas), 
                      alpha_generic_tilde = colMeans(PrelimPosteriors$alpha_generic_tilde), 
                      alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), 
                      local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
                      c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
                      alpha_hat_eij_tilde = colMeans(PrelimPosteriors$alpha_hat_eij_tilde), 
                      local_shrinkage_eij = colMeans(PrelimPosteriors$local_shrinkage_eij),
                      alpha_intra_tilde = colMeans(PrelimPosteriors$alpha_intra_tilde))
InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# Run the final fit of the model
FinalFit <- stan(file = FinalStanPath, data = FinalDataVec, iter = 3000,
                 chains = 3, init = InitVals, control = list(adapt_delta = 0.9))
Posteriors <- extract(FinalFit)

# Examine the same diagnostic plots as before to check for problems with the 
#       final model fit
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = "lambdas")
traceplot(FinalFit, pars = "alpha_generic")
traceplot(FinalFit, pars = "alpha_intra")
which(Inclusion_ij == 1)
traceplot(FinalFit, pars = "alpha_hat_ij")
which(Inclusion_eij == 1)
traceplot(FinalFit, pars = "alpha_hat_eij")

# Double check the autocorrelation
acf(Posteriors$lambdas[,1])
acf(Posteriors$lambdas[,2])
acf(Posteriors$alpha_generic[,1])
acf(Posteriors$alpha_generic[,2])
acf(Posteriors$alpha_intra[,1])
acf(Posteriors$alpha_intra[,2])
for(s in 1:S){
        if(Inclusion_ij[s] == 1){
                quartz()
                acf(Posteriors$alpha_hat_ij[,s])
        }
        if(Inclusion_eij[s] == 1){
                quartz()
                acf(Posteriors$alpha_hat_eij[,s])
        }
}

# If the fit looks good, safe the final output here
FitFileName <- paste(here("BH_simulations/Main/StanFits/"), FilePrefix, "FinalFit_R1.rdata", sep = "")
save(FinalFit, Posteriors, Inclusion_ij, Inclusion_eij, file = FitFileName)


########### Now calculate the values for the posterior predictive check and other 
#       assessments of parameter accuracy
PostLength <- length(Posteriors$alpha_generic[,1])
# calculate the posterior distributions of the interaction coefficients and lambdas
alpha_eij <- array(NA, dim = c(PostLength, N_ppc, S))
lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:N_ppc){
        lambda_ei[,i] <- exp(Posteriors$lambdas[,1] + Posteriors$lambdas[,2]*env_ppc[i])
        for(s in 1:S){
                if(s == Focal){
                        alpha_eij[,i,s] <- exp(Posteriors$alpha_intra[,1] + Posteriors$alpha_intra[,2] * env_ppc[i])
                }else{
                        alpha_eij[,i,s] <- exp(Posteriors$alpha_generic[1] + Inclusion_ij[s] * Posteriors$alpha_hat_ij[,s] +
                                (Posteriors$alpha_generic[,2] + Inclusion_eij[s] * Posteriors$alpha_hat_eij[,s]) * env_ppc[i])
                }
        }
}

# use the above quantities to calculate the posterior prediction intervals for the new data
Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
for(i in 1:PostLength){
     for(j in 1:N_ppc){
          SigmaTerm <- sum(alpha_eij[i,j,] * SpMatrix_ppc[j,])
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

# Now calculate the parameter deviations, starting with lambdas
LambdaEsts <- matrix(data = NA, nrow = 3, ncol = 2)
LambdaEsts[1,1] <- mean(Posteriors$lambdas[,1] - TrueVals$lambda.mean[Focal])
LambdaEsts[2:3,1] <- hdi(Posteriors$lambdas[,1] - TrueVals$lambda.mean[Focal])
LambdaEsts[1,2] <- mean(Posteriors$lambdas[,2] - TrueVals$lambda.env[Focal])
LambdaEsts[2:3,2] <- hdi(Posteriors$lambdas[,2] - TrueVals$lambda.env[Focal])

# Now calculate the alpha values
AlphaEsts <- array(data = NA, dim = c(3,2,S))
for(s in 1:S){
        if(s == Focal){
                # First the intercept
                Intercept <- Posteriors$alpha_intra[,1]
                AlphaEsts[1,1,s] <- mean(Intercept - TrueAlphaMeans[s])
                AlphaEsts[2:3,1,s] <- hdi(Intercept - TrueAlphaMeans[s])
                # Now the slope
                Slope <- Posteriors$alpha_intra[,2]
                AlphaEsts[1,2,s] <- mean(Slope - TrueAlphaSlopes[s])
                AlphaEsts[2:3,2,s] <- hdi(Slope - TrueAlphaSlopes[s])
        }else{
                # First the intercept
                Intercept <- Posteriors$alpha_generic[,1] + Posteriors$alpha_hat_ij[,s] * Inclusion_ij[s]
                AlphaEsts[1,1,s] <- mean(Intercept - TrueAlphaMeans[s])
                AlphaEsts[2:3,1,s] <- hdi(Intercept - TrueAlphaMeans[s])
                # Now the slope
                Slope <- Posteriors$alpha_generic[,2] + Posteriors$alpha_hat_eij[,s] * Inclusion_eij[s]
                AlphaEsts[1,2,s] <- mean(Slope - TrueAlphaSlopes[s])
                AlphaEsts[2:3,2,s] <- hdi(Slope - TrueAlphaSlopes[s])
        }
}

# Finally, calculate the "true" generic alpha that the model is attempting to estimate
GenericIntercepts <- 1 - Inclusion_ij
GenericSlopes <- 1 - Inclusion_eij
GenericIntercepts[Focal] <- 0
GenericSlopes[Focal] <- 0
TrueGenericIntercept <- log(sum(colSums(SpMatrix) * GenericIntercepts * exp(TrueAlphaMeans)) / sum(colSums(SpMatrix) * GenericIntercepts))
Totals <- colSums(SpMatrix) * GenericSlopes
TrueGenericSlope <- mean(TrueAlphaSlopes*GenericSlopes*(Totals/max(Totals)))

# Finally, save all the necessary results for the figures
FileName <- paste(here("BH_simulations/Main/StanFits/"), FilePrefix, "GraphStuff_R1.rdata", sep = "")
save(PredVals, Growth_ppc, LambdaEsts, AlphaEsts, Inclusion_eij, Inclusion_ij,
     TrueGenericIntercept, TrueGenericSlope,
     file = FileName)

