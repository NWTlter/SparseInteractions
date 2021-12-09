# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript


rm(list = ls())
library(tidyverse)
library(rstan)
library(here)

rm(list = ls(all.names = TRUE))

# ##  Print lots of things (TRUE), or not (FALSE)?
# verbose <- FALSE
# 
# source("utility_functions/utility_functions_all.R")
# 
# 
# ####
# ####  PACKAGES ----
# ####
# library(tidyverse) # data munging
# library(dplyr)     # more data munging
# library(codyn)     # community dynamics metrics
# library(vegan)     # diversity metrics
# library (ggpmisc) # to get the pvalues and slopes to print on the figs
# 
# ####
# ####  FILE NAMES, PATHS, READ DATA, ETC. ----
# ####
# #data_path       <- "../data/"
# #raw_cover_fname <- "filename.csv"
# #raw_anpp_fname  <- "filename.csv"
# #species_cover   <- read.csv(paste0(data_path, raw_cover_fname)) %>% select(-block)
# species_cover <- getTabular(93)
# 
# #but this is fast..need to look up the ID for pp
# anpp <- read_csv("http://niwot.colorado.edu/data_csvs/saddgrid_npp.hh.data.csv",
#                  guess_max = 10000
# )
# 
# 
# #Data Cleaning####
# #anpp
# #unique(anpp$year)
# 
# # use only 'A' or 'none' subsample as most consistent across time; remove snow
# # fence plots and shrub tundra for climate analyses
# anpp <-
#         anpp %>%
#         filter(subsample == "none" | subsample == "A") #%>%
# #filter(veg_class != "SF") %>%
# #filter(veg_class != "ST")
# 
# 
# #subset to veg classes
# veg_class <-
#         anpp %>%
#         filter(year == 2010) %>%
#         select(grid_pt, veg_class) %>%
#         distinct() %>%
#         rename (plot = grid_pt)
# 
# ####
# ####  1. AVERAGED SPECIES SYNCHRONY ----
# ####
# 
# #cleaning prep script copying M. Oldfather's NWT_LTER_Saddle_composition.R
# # read in community data file
# 
# # look at data columns
# 
# # 2RF = rock fragments, 2LICHN = lichen, 2X = registration marker, 2LTR = Litter,  2BARE = bare ground, 2FORB = unknown forb, 2GRAM = unknown graminoid, 2MOSS = moss etc ... will need ti remove the unknown species as well as seperate out lichen, moss, rock, litter (everything starting with 2 is a non-species)
# non_plant <- c("2RF", "2LICHN", "2X","2LTR","2BARE", "2HOLE", "2MOSS", "2SCATE")
# unknown_species <- c("2FORB","2MOSS","2GRAM", "2UNKSC", "POA","2UNK", "CAREX", "2COMP") 
# 
# # keep only bottom and top hits to allow for consistency across time
# species_cover <-
#         species_cover %>% 
#         filter(hit_type == "bottom" | hit_type == "top")
# 
# # change year = 1996 to 1995 for plot 37 (typo in dataset)
# species_cover[species_cover$year == 1996, "year"] <- 1995
# 
# # Additional subsetting --> top-most hit at each hit in each year (e.g. bottom hit if there is only a bottom hit, otherwise, top hit)
# # checking...
# 
# species_cover <- species_cover %>% 
#         arrange(year, plot, x, y, desc(hit_type)) %>% 
#         group_by(year, plot, x,y) %>% 
#         slice(1)
# 
# # aggregate data_file to calculate the number of hits (adds top and bottom together)
# species_cover <- species_cover %>% 
#         group_by(year, plot, USDA_code, USDA_name) %>% 
#         summarise(hits = n()) %>%
#         left_join(., veg_classes)
# 
# #most are ST or SF (very poorly represented)
# # but also a fair bit of barren
# # View (species_cover %>%
# #         ungroup() %>%
# #         filter(is.na(veg_class)) %>%
# #   select(plot) %>%
# #   distinct())
# 
# #73 isn't actually fully bare
# # maybe was supposed to be FF or??
# # but anyhow, drop for now
# species_cover = species_cover %>%
#         mutate(veg_class = ifelse (plot %in% c(1,41,42, 13, 71, 73, 301), 'B', veg_class)) %>%
#         filter(!is.na(veg_class)) %>%
#         filter(!veg_class %in% c('B', 'SF', 'ST'))
# 
# 
# species_cover <- species_cover %>%
#         ungroup()%>%
#         filter(!grepl('^2', USDA_code)) %>% #drop unks and not-plants
#         filter(!USDA_code %in% unknown_species) %>% #drop unks and not-plantsunknown_species
#         rename(plot_id = plot,
#                calendar_year = year,
#                cover = hits,
#                species = USDA_name) %>%
#         select(-USDA_code) %>%
#         pivot_wider(., names_from = species, values_from = cover,
#                     values_fill = 0)
# 
# write.csv(species_cover, file = "Empirical/saddle_grid.csv")

saddle_temp <- read_csv("Empirical/crall21xhcn_ctw_2020.csv")

#calculate GDD per year
gdd <- saddle_temp %>% 
        group_by(Year) %>% 
        mutate(TMEAN = (TMAX + TMIN)/2) %>% 
        mutate(mean_t = mean(TMEAN[Month %in% c(5:9)])) %>% 
        mutate(max_t = max(TMEAN)) %>% 
        filter(TMEAN > 5) %>% 
        mutate(gdd = sum(TMEAN)) %>% 
        filter(Month < 6) %>% 
        mutate(sp_gdd = sum(TMEAN)) %>% 
        select(Year, gdd, sp_gdd, mean_t, max_t) %>% 
        distinct() %>% 
        rename("year" = Year) %>% 
        filter(year > 2009)



#load in data #
saddle <- read.csv(here("Empirical/Saddle_grid.csv")) %>% 
        filter(veg_class %in% c("MM", "WM")) 

#year = 2017
#year1 = year + 1
focal_sp = "Geum rossii"
mod_type = "moist_cool"
FocalPrefix <- "GEUROS"
 
#low temp years
years = c(2011,2014,2016,2017, 2019)
#high temp years
#years = c(2010,2012,2013,2015,2018)

sad1 <- saddle %>% 
        filter(calendar_year %in% years)
     
sad2 <- saddle %>% 
     filter(calendar_year %in% (years+1)) %>%
     mutate(calendar_year = calendar_year - 1) %>% 
     select(Deschampsia.cespitosa, plot_id, calendar_year) %>% 
     rename("t1" = Deschampsia.cespitosa) 

abund <- sad1 %>% 
     left_join(sad2) %>% 
     mutate(t1 = ifelse(is.na(t1), 0, t1)) #%>% 
     #    select(-"NA")

# env <- read.csv(here("Empirical/snowmelt_est.csv")) %>% 
#      select("PLOT" = grid_pt, Estimate, "YEAR" = year) %>%
#         filter(YEAR > 2009) %>% 
#         mutate(YEAR = YEAR - 1)
        

all_dat <- abund #%>% 
     # left_join(env) 


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#FocalLetter <- "W" # "W" or "A"
 # "Waitzia" or "ARCA"
#FocalSpecies <- "Waitzia.acuminata" # "Waitzia.acuminata" or "Arctotheca.calendula"

#EnvCov <- "Estimate" # "Phos" or "Shade"
#EnvCol <- 51  # 72 for Canopy or 71 for Phosphorous

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
Nt <- as.integer(all_dat$Geum.rossii.var..turbinatum)
#reserve <- as.integer(as.factor(SpDataFocal$Reserve.x))
#reserve <- as.integer(as.factor(all_dat$Community_type))
#env <- as.vector(scale(SpDataFocal[,EnvCol]))
# env <- as.vector(scale(all_dat[,EnvCol]))

# Now calculate the total number of species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpAbunds <- all_dat[,5:110]
AllSpNames <- names(all_dat[,5:110])
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
Intra <- ifelse(SpNames == "Geum.rossii.var..turbinatum", 1, 0)

# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4
# DataVec <- c("N", "S", "Nt", "Fecundity", "reserve", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")


DataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "Intra", "tau0", "slab_scale", "slab_df")
Ntp1 <- Fecundity

# Now run a perliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = here("BH_simulations/Main/StanCode/Prelim_monoLambda_envAlpha_noenv.stan"), data = DataVec, iter = 1000, chains = 3, control = list(adapt_delta = 0.99))

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

DataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "Intra", "Inclusion_ij", "Inclusion_eij")

FinalFit <- stan(file = here("BH_simulations/Main/StanCode/Final_monoLambda_envAlpha_noenv.stan"), data = DataVec, iter = 3000, chains = 3, control = list(adapt_delta = 0.9))
FinalPosteriors <- rstan::extract(FinalFit)

# Diagnostic figures
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra"))
# acf(FinalPosteriors$lambdas[,1,1])
# acf(FinalPosteriors$lambdas[,1,2])
# acf(FinalPosteriors$lambdas[,2,1])
# acf(FinalPosteriors$lambdas[,2,2])
# acf(FinalPosteriors$alpha_generic[,1])
# acf(FinalPosteriors$alpha_generic[,2])
# acf(FinalPosteriors$alpha_intra[,1])
# acf(FinalPosteriors$alpha_intra[,2])

FileName <- paste(here("Empirical/StanFits/"), "/", FocalPrefix, "_", mod_type, "_FinalFit.rdata", sep = "")
save(FinalFit, SpNames, N, S, Fecundity, SpMatrix, Inclusion_ij,
     Inclusion_eij, tau0, slab_scale, slab_df, Intra, file = FileName)




#### output figure ####
library(bayestestR)
out <- data.frame(species = NULL, mod_tyoe = NULL, lambda_mean = NULL, lambda_ci = NULL, alpha.gen_mean = NULL, alpha.gen_ci = NULL, alpha.intra_mean = NULL, alpha.intra_ci = NULL)

for(j in c("KOBMYO", "CAMROT", "ARCSCO", "DESCES", "CHIJAM", "GEUROS")){
for(i in c("moist_cool", "moist_warm", "dry_warm", "dry_cool")){
        FileName <- paste(here("Empirical/StanFits/"), "/", j, "_", i, "_FinalFit.rdata", sep = "")
        
        load(FileName)
        FinalPosteriors <- rstan::extract(FinalFit)
        
        dat <- data.frame(species = j, 
                          mod_type = i, 
                          lambda_mean = mean(FinalPosteriors$lambdas), 
                          lambda_ci = ci(FinalPosteriors$lambdas), 
                          alpha.gen_mean = mean(exp(FinalPosteriors$alpha_generic)), 
                          alpha.gen_ci = ci(exp(FinalPosteriors$alpha_generic)), 
                          alpha.intra_mean = mean(exp(FinalPosteriors$alpha_intra)), 
                          alpha.intra_ci = ci(exp(FinalPosteriors$alpha_intra)))
        
        out <- bind_rows(out, dat)
        
}}

out %>% 
        pivot_longer(-c(species, mod_type), values_to = "values", names_to = "names") %>% 
        separate(names, sep = "_", into = c("variable", "type", "type2")) %>% 
        mutate(type = ifelse(!is.na(type2), type2, type)) %>% 
        select(-type2) %>% 
        pivot_wider(names_from = "type", values_from = "values") %>%
        separate(mod_type, sep = "_", into  = c("moisture", "temp")) %>% 
        mutate(variable = ifelse(variable == "alpha.gen", "alpha.inter", variable)) %>% 
        mutate(variable = factor(variable, levels = c("lambda", "alpha.inter", "alpha.intra"))) %>% 
        mutate(temp = factor(temp, levels = c("warm", "cool"))) %>%
        filter(moisture != "moist" | species != "CAMROT") %>% 
        filter(moisture != "dry" | species != "CHIJAM") %>% 
        ggplot(aes(x = moisture, y = mean, color = temp)) +
        geom_point(position = position_dodge(0.5), size = 2) +
        geom_errorbar(aes(x = moisture, ymin = low, ymax = high, color = temp), position = position_dodge(0.5), width = 0.1, size = 1) +
        facet_grid(variable~species, scales = "free_y") +
        theme(text = element_text(size = 15))

## CAMROT is not present in moist plots


out %>% 
        pivot_longer(-c(species, mod_type), values_to = "values", names_to = "names") %>% 
        separate(names, sep = "_", into = c("variable", "type", "type2")) %>% 
        mutate(type = ifelse(!is.na(type2), type2, type)) %>% 
        select(-type2) %>% 
        pivot_wider(names_from = "type", values_from = "values") %>%
        separate(mod_type, sep = "_", into  = c("moisture", "temp")) %>% 
        mutate(variable = ifelse(variable == "alpha.gen", "alpha.inter", variable)) %>% 
        mutate(variable = factor(variable, levels = c("lambda", "alpha.inter", "alpha.intra"))) %>% 
        mutate(temp = factor(temp, levels = c("warm", "cool"))) %>% 
        filter(moisture != "moist" | species != "CAMROT") %>% 
        filter(moisture != "dry" | species != "CHIJAM") %>% 
        filter(moisture != "dry" | species != "ARCSCO") %>% 
        ggplot(aes(x = moisture, y = mean, color = temp)) +
        geom_point(position = position_dodge(0.5), size = 2) +
        geom_errorbar(aes(x = moisture, ymin = low, ymax = high, color = temp), position = position_dodge(0.5), width = 0.1, size = 1) +
        facet_grid(variable~species, scales = "free_y") +
        theme(text = element_text(size = 15))


out %>% 
        mutate(intra_inter_ratio = alpha.intra_mean/alpha.gen_mean) %>% 
        select(species, mod_type, intra_inter_ratio) %>% 
        separate(mod_type, sep = "_", into = c("moisture", "temp")) %>%         filter(moisture != "moist" | species != "CAMROT") %>% 
        filter(moisture != "dry" | species != "CHIJAM") %>% 
        filter(moisture != "dry" | species != "ARCSCO") %>% 
        mutate(temp = factor(temp, levels = c("warm", "cool"))) %>%
        ggplot(aes(x = moisture, y = intra_inter_ratio, color = temp)) +
        geom_point(size = 3) +
        geom_hline(aes(yintercept = 1)) +
        facet_wrap(~species)+
        theme(text = element_text(size = 15))

