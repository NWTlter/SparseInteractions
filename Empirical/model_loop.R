# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

rm(list = ls())
library(tidyverse)
library(rstan)
library(here)

#community data
saddle <- read.csv(here("Empirical/Saddle_grid_clean.csv")) %>% 
        select(SPP, Raw_Abund, PLOT, YEAR) %>% 
        distinct() 

spp <- saddle %>% #pick species over 300 hits
     group_by(SPP) %>% 
     summarize(n = n()) %>% 
     filter(n > 300)

#set prior characteristics
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4

#other settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



for(i in 2010){
     for(j in spp$SPP){
          year = i
          year1 = year + 1
          focal_sp = j
          
          sad1 <- saddle %>% 
               filter(YEAR == year) %>% 
               select(-YEAR) %>% 
               pivot_wider(names_from = SPP, 
                           values_from = Raw_Abund, 
                           values_fill = 0)
          
          sad2 <- saddle %>% 
               filter(YEAR == year1) %>%
               filter(SPP == focal_sp) %>% 
               select(Raw_Abund, PLOT) %>% 
               rename("t1" = Raw_Abund) 
          
          abund <- sad1 %>% 
               left_join(sad2) %>% 
               mutate(t1 = ifelse(is.na(t1), 0, t1)) 
          
          env <- read.csv(here("Empirical/snowmelt_est.csv")) %>% 
               filter(year == year1) %>% 
               select("PLOT" = grid_pt, Estimate)
          
          all_dat <- abund %>% 
               left_join(env) 

          EnvCov <- "Estimate"
          EnvCol <- as.integer(ncol(all_dat))    
     
          N <- as.integer(nrow(all_dat))
          Ntp1 <- as.integer(all_dat$t1)
          
          Nt <- all_dat %>% select(j)
          Nt <- as.integer(unlist(Nt))
          
          env <- as.vector(scale(all_dat[,EnvCol]))
     
          AllSpAbunds <- all_dat[,2:(ncol(all_dat)-2)]
          AllSpNames <- names(all_dat[,2:(ncol(all_dat)-2)])
          SpTotals <- colSums(AllSpAbunds)
          SpToKeep <- SpTotals > 0
          S <- sum(SpToKeep)
         
          SpMatrix <- AllSpAbunds[SpToKeep]
          
          SpNames <- AllSpNames[SpToKeep]
          Intra <- ifelse(SpNames == focal_sp, 1, 0)
          
          DataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "tau0", "slab_scale", "slab_df")
          
          PrelimFit <- stan(file = here("BH_simulations/Main/StanCode/Prelim_monoLambda_envAlpha_quad.stan"), data = DataVec, iter = 3000, chains = 3, control = list(adapt_delta = 0.99))
          PrelimPosteriors <- rstan::extract(PrelimFit)
          
          Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = S)
          Inclusion_eij <- matrix(data = 0, nrow = 1, ncol = S)
          IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade

     for(s in 1:S){
                    Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
                    Ints_eij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = IntLevel)
                    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
                         Inclusion_ij[,s] <- 1
                    }
               
          }
          sum(Inclusion_ij)
          sum(Inclusion_eij)
          
          Inclusion_ij <- as.vector(Inclusion_ij)
          Inclusion_eij <- as.vector(Inclusion_eij)
          
          DataVec <- c("N", "S", "Nt", "Ntp1", "SpMatrix", "env", "Intra", "Inclusion_ij", "Inclusion_eij")
          
          FinalFit <- stan(file = here("BH_simulations/Main/StanCode/Final_monoLambda_envAlpha_quad.stan"), data = DataVec, iter = 3000, chains = 3, control = list(adapt_delta = 0.99))
          FinalPosteriors <- rstan::extract(FinalFit)
          
          FileName <- paste(here("Empirical/StanFits/"),"/", EnvCov, "_", i, "_", j, "_FinalFit_quad.rdata", sep = "")
          save(FinalFit, PrelimFit, SpNames, N, S, Nt, Ntp1, SpMatrix, env, Inclusion_ij, Inclusion_eij, tau0, slab_scale, slab_df, Intra, file = FileName)

          FigName <- paste(here("Empirical/Figures/"),"/", EnvCov, "_", i, "_", j, "_FinalFit_quad.pdf", sep = "")  
          
          ObsSnow <- as.vector(scale(all_dat$Estimate))
          EnvLength <- 1000 
          PlotSnow <- seq(min(ObsSnow, na.rm = TRUE), max(ObsSnow, na.rm = TRUE), length.out = EnvLength)
          
          Post <- rstan::extract(FinalFit)
          Snow <- list(Post = Post, SpNames = SpNames, N = N, S = S,  SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij, Inclusion_eij = Inclusion_eij, Intra = Intra)
          #rm(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij, Inclusion_eij, tau0, slab_scale, slab_df, Intra, Post)
          
          LambdaPlotVals <- array(NA, dim = c(1, 1, EnvLength, 3)) 
          
          for(k in 1:EnvLength){
               for(l in 1){
                    SnowPost <- exp(Snow$Post$lambdas[,1] + Snow$Post$lambdas[,2] * PlotSnow[k] + Snow$Post$lambdas[,3] * PlotSnow[k]^2)
                    #ShadePost <- exp(Shade$Post$lambdas[,j,1] + Shade$Post$lambdas[,j,2] * PlotShade[i])
                    LambdaPlotVals[1,l,k,1] <- mean(SnowPost)
                    LambdaPlotVals[1,l,k,2:3] <- HDInterval::hdi(SnowPost)
                    #LambdaPlotVals[2,j,i,1] <- mean(ShadePost)
                    #LambdaPlotVals[2,j,i,2:3] <- HDInterval::hdi(ShadePost)
               }
          }
          #rm(SnowPost)
          
          GenericPlot <- array(NA, dim = c(1,3,EnvLength))
          IntraPlot <- array(NA, dim = c(1,3,EnvLength))
          SpecificPlot <- array(NA, dim = c(1,2, 3, EnvLength))
          
          SnowLegend <- c("Generic", "Intraspecific")
          
          for(m in 1:EnvLength){
               GenericPost <- exp(Snow$Post$alpha_generic[,1] + Snow$Post$alpha_generic[,2] * PlotSnow[m] + Snow$Post$alpha_generic[,3] * PlotSnow[m]^2)
               IntraPost <- exp(Snow$Post$alpha_intra[,1] + Snow$Post$alpha_intra[,2] * PlotSnow[m] + Snow$Post$alpha_intra[,3] * PlotSnow[m]^2)
              
               GenericPlot[1,1,m] <- mean(GenericPost)
               GenericPlot[1,2:3,m] <- HDInterval::hdi(GenericPost)
               IntraPlot[1,1,m] <- mean(IntraPost)
               IntraPlot[1,2:3,m] <- HDInterval::hdi(IntraPost)
          }
          library(RColorBrewer)
          AlphaRange <- c(0, max(GenericPlot[1,1,]))
          LambdaRange <- range(LambdaPlotVals)
          SnowRange <- range(ObsSnow, na.rm = TRUE)
          # ShadeRange <- range(ObsShade, na.rm = TRUE)
          DarkCols <- brewer.pal(n = 8, name = "Dark2")
          LambdaCols <- DarkCols[3:4]
          IntraCol <- DarkCols[5]
          AlphaCols <- DarkCols[6:7]
          LambdaLab <- expression(lambda["e,i"])
          AlphaLab <- expression(alpha["e,i,j"])
          
          pdf(file = FigName, width = 4, height = 6, onefile = FALSE, paper = "special")
          par(mfrow = c(2,1), oma = c(4,4,2,2), mar = c(2,2,2,2))
          # First, plot the lambda results for phosphorous
          plot(NA, NA, xlim = SnowRange, ylim = LambdaRange, main = "", xlab = "", ylab = LambdaLab,
               las = 1, xpd = NA, cex.lab = 1.5)
          axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
          axis(side = 2, at = seq(0, 20, by = 1), labels = FALSE, tcl = -0.25)
          for(n in 1){
               lines(x = PlotSnow, y = LambdaPlotVals[1,n,,1], col = LambdaCols[n], lwd = 1.5)
               lines(x = PlotSnow, y = LambdaPlotVals[1,n,,2], col = LambdaCols[n], lty = 2)
               lines(x = PlotSnow, y = LambdaPlotVals[1,n,,3], col = LambdaCols[n], lty = 2)
          }

          text(x = 0.95*SnowRange[1], y = 0.95*LambdaRange[2], 
               labels = expression(bold("a")), cex = 1.5)
     
          plot(NA, NA, xlim = SnowRange, ylim = AlphaRange, main = "", xlab = "Standardized Snow", 
               ylab = AlphaLab, las = 1, xpd = NA, cex.lab = 1.5)
          axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
          axis(side = 2, at = seq(0, 0.06, by = 0.005), labels = FALSE, tcl = -0.25)
          # Add the generic polygon and line
          xCoords <- c(PlotSnow, PlotSnow[EnvLength:1])
          yCoords <- c(GenericPlot[1,2,], GenericPlot[1,3,EnvLength:1])
          polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
          lines(x = PlotSnow, y = GenericPlot[1,1,], lty = 1, lwd = 1.5)
          # Add the lines for the intraspecific term
          lines(x = PlotSnow, y = IntraPlot[1,1,], lwd = 1.5, col = IntraCol)
          #lines(x = PlotSnow, y = IntraPlot[1,2,], lty = 2, col = IntraCol)
          #lines(x = PlotSnow, y = IntraPlot[1,3,], lty = 2, col = IntraCol)
          # Now the deviating alphas
          legend("topright", legend = SnowLegend, lty = 1, lwd = 1.5, col = c("black", IntraCol, AlphaCols),
                 bty = "n")
          text(x = 0.95*SnowRange[1], y = 0.95*AlphaRange[2], 
               labels = expression(bold("c")), cex = 1.5)
          dev.off()
          
     }
}









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
