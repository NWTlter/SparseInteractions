# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript


rm(list = ls())
library(tidyverse)
library(rstan)
library(here)

#load in data #
saddle <- read.csv(here("Empirical/Saddle_grid_clean.csv")) %>% 
     select(SPP, Raw_Abund, PLOT, YEAR) %>% 
     distinct() 

year = 2017
year1 = year + 1
focal_sp = "DESCES"
 
sad1 <- saddle %>% 
     filter(YEAR == year) %>% 
     select(-YEAR) %>% 
     pivot_wider(names_from = SPP, values_from = Raw_Abund, values_fill = 0)
     
sad2 <- saddle %>% 
     filter(YEAR == year1) %>%
     filter(SPP == focal_sp) %>% 
     select(Raw_Abund, PLOT) %>% 
     rename("t1" = Raw_Abund) 

abund <- sad1 %>% 
     left_join(sad2) %>% 
     mutate(t1 = ifelse(is.na(t1), 0, t1)) %>% 
        select(-"NA")

env <- read.csv(here("Empirical/snowmelt_est.csv")) %>% 
     filter(year == year1) %>% 
     select("PLOT" = grid_pt, Estimate)

all_dat <- abund %>% 
     left_join(env) 


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#FocalLetter <- "W" # "W" or "A"
FocalPrefix <- "DESCES" # "Waitzia" or "ARCA"
#FocalSpecies <- "Waitzia.acuminata" # "Waitzia.acuminata" or "Arctotheca.calendula"

EnvCov <- "Estimate" # "Phos" or "Shade"
EnvCol <- 59  # 72 for Canopy or 71 for Phosphorous

# Load in the data and subset out the current focal species.
#SpData <- read.csv(here("Empirical/water_full_env.csv"))
#SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
#SpData <- na.omit(SpData) 
#FocalLetter
#SpDataFocal <- subset(SpData, Focal.sp.x == FocalLetter)



# Next continue to extract the data needed to run the model. 
#N <- as.integer(nrow(SpDataFocal))
N <- as.integer(nrow(all_dat))
#Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
Fecundity <- as.integer(all_dat$t1)
Nt <- as.integer(all_dat$GEUROS)
#reserve <- as.integer(as.factor(SpDataFocal$Reserve.x))
#reserve <- as.integer(as.factor(all_dat$Community_type))
#env <- as.vector(scale(SpDataFocal[,EnvCol]))
env <- as.vector(scale(all_dat[,EnvCol]))

# Now calculate the total number of species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpAbunds <- all_dat[,2:57]
AllSpNames <- names(all_dat[,2:57])
SpTotals <- colSums(AllSpAbunds)
SpToKeep <- SpTotals > 0
S <- sum(SpToKeep)
# SpMatrix <- matrix(NA, nrow = N, ncol = S)
# i <- 1
# for(s in 1:ncol(AllSpAbunds)){
#      if(SpToKeep[s] == T){
#           SpMatrix[,i] <- AllSpAbunds[,s]
#           i <- i + 1
#      }
# }

SpMatrix <- AllSpAbunds[SpToKeep]

SpNames <- AllSpNames[SpToKeep]
Intra <- ifelse(SpNames == focal_sp, 1, 0)

# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4
# DataVec <- c("N", "S", "Nt", "Fecundity", "reserve", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")

DataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
Ntp1 <- Fecundity

# Now run a perliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = here("BH_simulations/Main/StanCode/Prelim_monoLambda_envAlpha_quad.stan"), data = DataVec, iter = 3000, chains = 3, control = list(adapt_delta = 0.99))

# PrelimFit <- stan(file = here("Empirical/StanCode/BH_FH_Preliminary.stan"), data = DataVec, iter = 3000, 
#                   chains = 1)
PrelimPosteriors <- rstan::extract(PrelimFit)

##### Diagnostic plots
# First check the distribution of Rhats and effective sample sizes
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
# Next check the correlation among key model parameters and identify any
#       divergent transitions
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
# Finally, check for autocorrelation in the posteriors of key model parameters
# acf(PrelimPosteriors$lambdas[,1,1])
# acf(PrelimPosteriors$lambdas[,1,2])
# acf(PrelimPosteriors$lambdas[,2,1])
# acf(PrelimPosteriors$lambdas[,2,2])
# acf(PrelimPosteriors$alpha_generic[,1])
# acf(PrelimPosteriors$alpha_generic[,2])
# acf(PrelimPosteriors$alpha_intra[,1])
# acf(PrelimPosteriors$alpha_intra[,2])

#### If the diagnostic plots don't reveal any problems wiht the model fit, now
#       move on to determining which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = S)
Inclusion_eij <- matrix(data = 0, nrow = 1, ncol = S)
IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
for(i in 1){
     for(s in 1:S){
          Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
          Ints_eij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = IntLevel)
          if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
               Inclusion_ij[i,s] <- 1
          }
     }
}
sum(Inclusion_ij)
sum(Inclusion_eij)

Inclusion_ij <- as.vector(Inclusion_ij)
Inclusion_eij <- as.vector(Inclusion_eij)

DataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij", "Inclusion_eij")

FinalFit <- stan(file = here("BH_simulations/Main/StanCode/Final_monoLambda_envAlpha_quad.stan"), data = DataVec, iter = 3000, chains = 3, control = list(adapt_delta = 0.99))
FinalPosteriors <- rstan::extract(FinalFit)

# Diagnostic figures
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
acf(FinalPosteriors$lambdas[,1,1])
acf(FinalPosteriors$lambdas[,1,2])
acf(FinalPosteriors$lambdas[,2,1])
acf(FinalPosteriors$lambdas[,2,2])
acf(FinalPosteriors$alpha_generic[,1])
acf(FinalPosteriors$alpha_generic[,2])
acf(FinalPosteriors$alpha_intra[,1])
acf(FinalPosteriors$alpha_intra[,2])

FileName <- paste(here("Empirical/StanFits/"), FocalPrefix, "_", EnvCov, "_FinalFit_quad.rdata", sep = "")
save(FinalFit, SpNames, N, S, Fecundity, SpMatrix, env, Inclusion_ij,
     Inclusion_eij, tau0, slab_scale, slab_df, Intra, file = FileName)

