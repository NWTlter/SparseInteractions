# File paths are assuming that the wd is where this source file is.
# These paths are within the "Box" folder (this file's location), and may need to be updated
# to the user's file structure

rm(list = ls())
library(here)
library(tidyverse)
library(patchwork)
library(HDInterval)

# plot theme
theme_cw <- function () { 
     theme_bw(base_size=12) %+replace% 
          theme(
               panel.background = element_blank(), 
               plot.background = element_blank(), 
               axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
               panel.grid.minor = element_blank(), 
               panel.grid.major = element_blank(),
               legend.background = element_blank(), 
               legend.key = element_blank(),
               strip.background = element_blank(), 
               complete = TRUE
          )
}

### Box Figure part 1: different formulas for lambda ---------------

# monotonic lambda data
# Model fits generated by monoLambda_constAlpha_modelFit.R file
# These paths are within the "Box" folder (this file's location), and may need to be updated
# to the user's file structure
load(here("BH_simulations/Box/StanFits/monoLambda_constAlpha/N50_FinalFit.rdata"))
PostLength <- length(FinalPosteriors$alpha_generic)
Focal <- 1

# individual draws from the posterior
n = 50
env.seq <- seq(-2, 2, length.out = n)
lambda.ei <- tibble(env.seq = rep(env.seq, times = PostLength), 
                    post = rep(1:PostLength, each = n), 
                    lambda.int = rep(FinalPosteriors$lambdas[,1], each = n),
                    lambda.slope = rep(FinalPosteriors$lambdas[,2], each = n))
lambda.ei$lambda.mono <- lambda.ei %>% with(exp(lambda.int + lambda.slope*env.seq))
lambda.ei$post <- as.factor(lambda.ei$post)

# sample 1000 from this to make graphing faster
sample.use <- sample(1:PostLength, size = 1000)
lambda.ei.small <- filter(lambda.ei, post %in% sample.use)
lambda.ei.small$ind <- 'ind'

# average from the posterior (could leave this out)
post.mean <- mean(FinalPosteriors$lambdas[,1])
post.env <- mean(FinalPosteriors$lambdas[,2])
lambda.post <- tibble(env.seq, post = as.factor(0), ind = 'post.mean', 
                      lambda.int = rep(post.mean, each = n),
                      lambda.slope = rep(post.env, each = n))
lambda.post$lambda.mono <- lambda.post %>% with(exp(lambda.int + lambda.slope*env.seq))


# true value
TrueVals <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/parameters_perturb2.csv"))
lambda.mean <- with(TrueVals, lambda.mean[species == Focal])
lambda.env <- with(TrueVals, lambda.env[species == Focal])
lambda.true <- tibble(env.seq, post = as.factor(0), ind = 'true', 
                      lambda.int = rep(lambda.mean, each = n),
                      lambda.slope = rep(lambda.env, each = n))
lambda.true$lambda.mono <- lambda.true %>% with(exp(lambda.int + lambda.slope*env.seq))

lambda.comp <- rbind(lambda.post, lambda.true)

# calculate comparisons
(lambda.mean - post.mean)/lambda.mean
(lambda.env - post.env)/lambda.env

p.mono <- ggplot(lambda.ei.small, aes(x = env.seq, y = lambda.mono)) + 
     geom_line(aes(group = post), alpha = 0.02, color = 'grey50') + 
     geom_line(data = lambda.comp, aes(color = ind, linetype = ind)) +
     scale_color_manual(values = c('black','red'), name = '', 
                        labels = c('Modeled','True')) + 
     scale_linetype_manual(values = c('dashed','solid'), name = '', 
                           labels = c('Modeled','True')) + 
     theme_cw() + 
     theme(legend.position = c(0.2, 0.9)) +
     ylab(expression(lambda[ei])) +
     xlab('Environment')

## Same process with optimum lambda
# Model fits generated by optLambda_constAlpha_modelFit.R file
# These paths are within the "Box" folder (this file's location), and may need to be updated
# to the user's file structure
load(here("BH_simulations/Box/StanFits/optLambda_constAlpha/N50_FinalFit.rdata"))
PostLength <- length(FinalPosteriors$alpha_generic)
Focal <- 8

# individual draws from the posterior
n = 50
env.seq <- seq(-2, 2, length.out = n)
lambda.ei <- tibble(env.seq = rep(env.seq, times = PostLength), 
                    post = rep(1:PostLength, each = n), 
                    lambda.max = rep(FinalPosteriors$lambda_max, each = n),
                    z.env = rep(FinalPosteriors$lambda_opt, each = n),
                    sigma.env = rep(FinalPosteriors$lambda_width, each = n))
lambda.ei$lambda.opt <- lambda.ei %>% 
     with(lambda.max * exp(-((z.env - env.seq)/(2*sigma.env))^2))
lambda.ei$post <- as.factor(lambda.ei$post)

# sample 1000 from this to make graphing faster
sample.use <- sample(1:PostLength, size = 1000)
lambda.ei.small <- filter(lambda.ei, post %in% sample.use)
lambda.ei.small$ind <- 'ind'

# average from the posterior (could leave this out)
post.max <- mean(FinalPosteriors$lambda_max)
post.z.env <- mean(FinalPosteriors$lambda_opt)
post.sigma.env <- mean(FinalPosteriors$lambda_width)
lambda.post <- tibble(env.seq, post = as.factor(0), ind = 'post.mean', 
                      lambda.max = rep(post.max, each = n),
                      z.env = rep(post.z.env, each = n),
                      sigma.env = rep(post.sigma.env, each = n))
lambda.post$lambda.opt <- lambda.post %>% 
     with(lambda.max * exp(-((z.env - env.seq)/(2*sigma.env))^2))


# true value
TrueVals <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/parameters_perturb_opt_const.csv"))
lambda.max <- with(TrueVals, lambda.max[species == Focal])
z.env <- with(TrueVals, z.env[species == Focal])
sigma.env <- with(TrueVals, sigma.env[species == Focal])
lambda.true <- tibble(env.seq, post = as.factor(0), ind = 'true', 
                      lambda.max = rep(lambda.max , each = n),
                      z.env = rep(z.env , each = n),
                      sigma.env = rep(sigma.env, each = n))
lambda.true$lambda.opt <- lambda.true %>% 
     with(lambda.max * exp(-((z.env - env.seq)/(2*sigma.env))^2))

# comparisons
(z.env - post.z.env)/z.env
(sigma.env - post.sigma.env)/sigma.env
(lambda.max - post.max)/lambda.max


p.opt <- ggplot(lambda.ei.small, aes(x = env.seq, y = lambda.opt)) + 
     geom_line(aes(group = post), alpha = 0.02, color = 'grey50') + 
     geom_line(data = lambda.true, color = 'red') +
     geom_line(data = lambda.post, color = 'black', linetype = 'dashed') + 
     theme_cw() + 
     ylab(expression(lambda[ei])) +
     xlab('Environment')

### Box figure part 2: different options for alpha------------
# Create lists for the results from different sizes of datasets
sample.size <- c(20, 50, 100, 200)
file.prefixes <- paste('N', sample.size, '_', sep = "")
n.samples <- length(file.prefixes)
alpha.type <- c('envAlpha', 'constAlpha')
Focal <- 1
S <- 15

alpha.fit <- tibble(sample.size = rep(sample.size, times = length(alpha.type)),
                    alpha.type = rep(alpha.type, each = n.samples),
                    n.ij = -1, n.eij = -1, n.total = -1,
                    dev.mean = -1, dev.sd = -1, dev.se = -1)

## env alpha fits
# Model fits used here generated by monoLambda_envAlpha_modelFit.R file
# These paths are nested within the "Box" folder (this file's location)
# and may need to be updated to the user's file structure

#FullSim <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/simulation_perturb2.csv"))
load('BH_simulations/test_multiple_simulations.rdata')
sim <- simulations[[8]]
FullSim <- sim$simulation

ppc_data <- subset(FullSim, (species == Focal) & (run > 200) & (time == 0) & (thinned == 0))
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

for(sn in 1:n.samples){
        ## number of samples included
        file.prelim <- paste(here("BH_simulations/Box/StanFits/monoLambda_envAlpha/"),
                             file.prefixes[sn],"PrelimFit.rdata", sep = "")
        load(file.prelim)
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
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'n.ij'] <- sum(Inclusion_ij)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'n.eij'] <- sum(Inclusion_eij)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'n.total'] <- sum(Inclusion_ij + Inclusion_eij)
        
        ## ppc deviations
        file.final <- paste(here("BH_simulations/Box/StanFits/monoLambda_envAlpha/"),
                            file.prefixes[sn],"FinalFit.rdata", sep = "")
        load(file.final)
        
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
        GrowthRMSE <- numeric(length = PostLength)
        for(i in 1:PostLength){
                for(j in 1:N_ppc){
                        SigmaTerm <- sum(alpha_eij[i,j,] * SpMatrix_ppc[j,])
                        Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
                        Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
                }
                deviation_sq <- (Growth_pred[i,] - Growth_ppc)^2
                GrowthRMSE[i] <- sqrt(sum(deviation_sq) / N_ppc)
        }
        
        # Calculate final fit results for the ppc
        PredVals <- colMeans(Growth_pred)
        pred.dev <- Growth_ppc - PredVals
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'dev.mean'] <- mean(abs(pred.dev))
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'dev.sd'] <- sd(abs(pred.dev))
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'dev.se'] <- sd(abs(pred.dev))/sqrt(length(pred.dev))
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'rmse'] <- mean(GrowthRMSE)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'rmse.low'] <- HDInterval::hdi(GrowthRMSE)[1]   
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'rmse.high'] <- HDInterval::hdi(GrowthRMSE)[2] 
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'tau'] <- mean(PrelimPosteriors$tau)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'tau.low'] <- HDInterval::hdi(PrelimPosteriors$tau)[1]   
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'envAlpha', 'tau.high'] <- HDInterval::hdi(PrelimPosteriors$tau)[2] 
}


## const alpha fits
# Model fits used here generated by monoLambda_constAlpha_modelFit.R file
# These paths are nested within the "Box" folder (this file's location)
# and may need to be updated to the user's file structure

#FullSim <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/simulation_perturb2_const.csv"))
FullSim <- read.csv(here("BH_simulations/Box/SimulationsDataFiles/simulation_const.csv")) %>%
        select(-X)

ppc_data <- subset(FullSim, (species == Focal) & (run > 200) & (time == 0) & (thinned == 0))
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
## number of samples included
for(sn in 1:n.samples){
        file.prelim <- paste(here("BH_simulations/Box/StanFits/monoLambda_constAlpha/"),
                             file.prefixes[sn],"PrelimFit.rdata", sep = "")
        load(file.prelim)
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
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'n.ij'] <- sum(Inclusion_ij)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'n.total'] <- sum(Inclusion_ij)
        
        ## ppc deviations
        file.final <- paste(here("BH_simulations/Box/StanFits/monoLambda_constAlpha/"),
                            file.prefixes[sn],"FinalFit.rdata", sep = "")
        load(file.final)
        
        PostLength <- length(FinalPosteriors$alpha_generic)
        # calculate the posterior distributions of the interaction coefficients and lambdas
        alpha_ij <- array(NA, dim = c(PostLength, N_ppc, S))
        lambda_ei <- matrix(NA, nrow = PostLength, ncol = N_ppc)
        for(i in 1:N_ppc){
                lambda_ei[,i] <- exp(FinalPosteriors$lambdas[,1] + FinalPosteriors$lambdas[,2]*env_ppc[i])
                for(s in 1:S){
                        if(s == Focal){
                                alpha_ij[,i,s] <- exp(FinalPosteriors$alpha_intra)
                        }else{
                                alpha_ij[,i,s] <- exp(FinalPosteriors$alpha_generic + Inclusion_ij[s] * FinalPosteriors$alpha_hat_ij[,s])
                        }
                }
        }
        
        # use the above quantities to calculate the posterior prediction intervals for the new data
        Growth_pred <- matrix(data = NA, nrow = PostLength, ncol = N_ppc)
        GrowthRMSE <- numeric(length = PostLength)
        
        for(i in 1:PostLength){
                for(j in 1:N_ppc){
                        SigmaTerm <- sum(alpha_ij[i,j,] * SpMatrix_ppc[j,])
                        Ntp1_pred <- Nt_ppc[j] * lambda_ei[i,j] / (1 + SigmaTerm)
                        Growth_pred[i,j] <- log((Ntp1_pred + 1)/Nt_ppc[j])
                }
                deviation_sq <- (Growth_pred[i,] - Growth_ppc)^2
                GrowthRMSE[i] <- sqrt(sum(deviation_sq) / N_ppc)
        }
        
        # Calculate final fit results for the ppc
        PredVals <- colMeans(Growth_pred)
        pred.dev <- Growth_ppc - PredVals
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'dev.mean'] <- mean(abs(pred.dev))
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'dev.sd'] <- sd(abs(pred.dev))
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'dev.se'] <- sd(abs(pred.dev))/sqrt(length(pred.dev))
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'rmse'] <- mean(GrowthRMSE)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'rmse.low'] <- HDInterval::hdi(GrowthRMSE)[1]   
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'rmse.high'] <- HDInterval::hdi(GrowthRMSE)[2]   
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'tau'] <- mean(PrelimPosteriors$tau)
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'tau.low'] <- HDInterval::hdi(PrelimPosteriors$tau)[1]   
        alpha.fit[alpha.fit$sample.size == sample.size[sn] & 
                          alpha.fit$alpha.type == 'constAlpha', 'tau.high'] <- HDInterval::hdi(PrelimPosteriors$tau)[2] 
}

## figure
# colors selected from the inauguration() package
InterceptCol <- "#f3c483"
SlopeCol <- "#5c1a33"

alpha.fit[alpha.fit == -1] <- NA
alpha.fit$sample.size.2 <- factor(alpha.fit$sample.size)
alpha.fit.long <- alpha.fit %>% 
        pivot_longer(cols = starts_with('n'), names_to = 'specific', values_to = 'number') %>%
        filter(!(alpha.type == 'constAlpha' & specific == 'n.total') & !(specific == 'n.eij'))

a.terms <- ggplot(alpha.fit.long, 
                  aes(x = sample.size, y = number, 
                      color = interaction(alpha.type, specific),
                      linetype = interaction(alpha.type, specific))) + 
        geom_line() +
        theme_cw() + 
        theme(legend.position = c(0.75, 0.2)) + 
        scale_color_manual(values = c(InterceptCol, SlopeCol, SlopeCol), name = "",
                           labels = c(expression(Simple~context:~alpha[ij]~pairs), 
                                      expression(Complex~context:~alpha[ij]~pairs),
                                      expression(Complex~context:~all~terms)),
                           guide = guide_legend(label.hjust = 0)) +
        scale_linetype_manual(values = c('solid','solid', 'dashed'), name = "",
                              labels = c(expression(Simple~context:~alpha[ij]~pairs), 
                                         expression(Complex~context:~alpha[ij]~pairs),
                                         expression(Complex~context:~all~terms)),
                              guide = guide_legend(label.hjust = 0)) +
        ylab('Number of non-generic terms') +
        xlab('Input data sample size')

a.post <- ggplot(alpha.fit, aes(x = sample.size, y = rmse, 
                                ymin = rmse.low, ymax = rmse.high,
                                color = alpha.type)) + 
        geom_point(position = position_dodge(5)) + 
        geom_errorbar(position = position_dodge(5), width = 0) +
        theme_cw() + 
        scale_color_manual(values = c(InterceptCol, SlopeCol), name = '',
                           labels = c('Simple context','Complex context'),
                           guide = guide_legend(label.hjust = 0)) +
        ylab('Root mean squared error') +
        xlab('Input data sample size') +
        theme(legend.position = c(0.8, 0.9))


## Not used: looking at how tau varies with sample size
a.tau <- ggplot(alpha.fit, aes(x = sample.size, y = tau, 
                               ymin = tau.low, ymax = tau.high,
                               color = alpha.type)) + 
        geom_point(position = position_dodge(5)) + 
        geom_errorbar(position = position_dodge(5), width = 0) +
        theme_cw() + 
        scale_color_manual(values = c(InterceptCol, SlopeCol), name = '',
                           labels = c('Simple context','Complex context'),
                           guide = guide_legend(label.hjust = 0)) +
        ylab('Tau') +
        xlab('Input data sample size') +
        theme(legend.position = c(0.8, 0.8))


### Box figure combining parts ---------------

patchwork <- (p.mono + p.opt) / (a.post + a.terms)

patchwork + 
        plot_annotation(tag_levels = 'a')  & 
        theme(plot.tag = element_text(face = "bold"))

# ggsave(filename = here('BH_simulations/Box/box_figure.pdf'), width = 10, height = 8, units = 'in')

