# This script will make the empirical figures for the York Gum Woodland data.
#    There will be two figures (1 per species), each with a plot of lambda and
#    alphas across both environmental gradients

rm(list = ls())
library(RColorBrewer)
library(here)

# Set the figure name
FigName <- here("Empirical/GEUROS2018_quad.pdf")

#load in data #
saddle <- read.csv(here("Empirical/Saddle_grid_clean.csv")) %>% 
     select(SPP, Raw_Abund, PLOT, YEAR) %>% 
     distinct() 

year = 2017
year1 = year + 1
focal_sp = "GEUROS"

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
     mutate(t1 = ifelse(is.na(t1), 0, t1))%>% 
     select(-"NA")

env <- read.csv(here("Empirical/snowmelt_est.csv")) %>% 
     filter(year == year1) %>% 
     select("PLOT" = grid_pt, Estimate)

all_dat <- abund %>% 
     left_join(env)



# Load in the empirical data and make objects for graphing
#FullData <- read.csv(here("Empirical/water_full_env.csv"))
ObsSnow <- as.vector(scale(all_dat$Estimate))
#ObsShade <- as.vector(scale(FullData$Canopy))
EnvLength <- 1000 
PlotSnow <- seq(min(ObsSnow, na.rm = TRUE), max(ObsSnow, na.rm = TRUE), length.out = EnvLength)
#PlotShade <- seq(min(ObsShade, na.rm = TRUE), max(ObsShade, na.rm = TRUE), length.out = EnvLength)

#ReserveNames <- c("DM", "MM")

# Load in the model fits for the current species
# Phosphorous
load(here("Empirical/StanFitsGEUROS_Estimate_FinalFit_quad.rdata"))
Post <- rstan::extract(FinalFit)
Snow <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij, Inclusion_eij = Inclusion_eij, Intra = Intra)
rm(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
   Inclusion_eij, tau0, slab_scale, slab_df, Intra, Post)
# Canopy cover
# load(here("Empirical/StanFits/WAAC_Shade_FinalFit.rdata"))
# Post <- rstan::extract(FinalFit)
# Shade <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
#               reserve = reserve, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij,
#               Inclusion_eij = Inclusion_eij, Intra = Intra)
# rm(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
#    Inclusion_eij, tau0, slab_scale, slab_df, Intra, Post)

# First, calculate the lambda values over both environmental variables
LambdaPlotVals <- array(NA, dim = c(1, 1, EnvLength, 3)) # Environmental co-variate, reserve, environmental sequence, metric (mean, lwr, upr)
for(i in 1:EnvLength){
     for(j in 1){
          SnowPost <- exp(Snow$Post$lambdas[,1] + Snow$Post$lambdas[,2] * PlotSnow[i] + Snow$Post$lambdas[,3] * PlotSnow[i]^2)
          #ShadePost <- exp(Shade$Post$lambdas[,j,1] + Shade$Post$lambdas[,j,2] * PlotShade[i])
          LambdaPlotVals[1,j,i,1] <- mean(SnowPost)
          LambdaPlotVals[1,j,i,2:3] <- HDInterval::hdi(SnowPost)
          #LambdaPlotVals[2,j,i,1] <- mean(ShadePost)
          #LambdaPlotVals[2,j,i,2:3] <- HDInterval::hdi(ShadePost)
     }
}
rm(SnowPost, ShadePost)

# Create objects to hold the mean and credible interval values for the intraspecific
#  and generic alpha values
GenericPlot <- array(NA, dim = c(1,3,EnvLength))
IntraPlot <- array(NA, dim = c(1,3,EnvLength))

# Create an object to hold the alpha values for any species identified by the sparse
#  modeling approach as deviating from the generic value
SpecificPlot <- array(NA, dim = c(1,2, 3, EnvLength))

# Create appropriate legends with the names for any deviating species
# Hyalosperma glutinosum and Schoenus nanus for phosphorous
# SnowLegend <- c("Generic", "Intraspecific",
#                 expression(paste(italic("IDK"), " IDK", sep = "")), 
#                 expression(paste(italic("IDK"), " IDK", sep = "")))

SnowLegend <- c("Generic", "Intraspecific")

# Hypochaeris glabra for shade
# ShadeLegend <- c("Generic", "Intraspecific",
#                  expression(paste(italic("H. glabra"), " (Bendering)", sep = "")))

# Now calculate the alpha_e,i,j values for intraspecific alphas (IntraPlot), any heterospecifics that
#  affect the focal species as an average heterospecific neighbor (GenericPlot), and for any species
#  that differ from that generic heterospecific effect (SpecificPlot)
for(i in 1:EnvLength){
     # First phosphorous
     # Calculate the posterior distribution for the relevant alpha_e,i,j at this particular environmental value
     GenericPost <- exp(Snow$Post$alpha_generic[,1] + Snow$Post$alpha_generic[,2] * PlotSnow[i] + Snow$Post$alpha_generic[,3] * PlotSnow[i]^2)
     #Sp1Post <- exp(Snow$Post$alpha_generic[,1] + Snow$Post$alpha_hat_ij[,2,17] + Snow$Post$alpha_generic[,2] * PlotSnow[i])
     #Sp2Post <- exp(Snow$Post$alpha_generic[,1] + Snow$Post$alpha_hat_ij[,1,37] + Snow$Post$alpha_generic[,2] * PlotSnow[i])
     IntraPost <- exp(Snow$Post$alpha_intra[,1] + Snow$Post$alpha_intra[,2] * PlotSnow[i] + Snow$Post$alpha_intra[,3] * PlotSnow[i]^2)
     # Now store the mean and 95% CI for these posteriors
     GenericPlot[1,1,i] <- mean(GenericPost)
     GenericPlot[1,2:3,i] <- HDInterval::hdi(GenericPost)
     IntraPlot[1,1,i] <- mean(IntraPost)
     IntraPlot[1,2:3,i] <- HDInterval::hdi(IntraPost)
     #SpecificPlot[1,1,1,i] <- mean(Sp1Post)
     #SpecificPlot[1,1,2:3,i] <- HDInterval::hdi(Sp1Post)
     #SpecificPlot[1,2,1,i] <- mean(Sp2Post)
     #SpecificPlot[1,2,2:3,i] <- HDInterval::hdi(Sp2Post)
     
     # Now shade, calculated in the same way as for phosphorous
     # GenericPost <- exp(Shade$Post$alpha_generic[,1] + Shade$Post$alpha_generic[,2] * PlotShade[i])
     # Sp1Post <- exp(Shade$Post$alpha_generic[,1] + Shade$Post$alpha_hat_ij[,1,19] + Shade$Post$alpha_generic[,2] * PlotShade[i])
     # IntraPost <- exp(Shade$Post$alpha_intra[,1] + Shade$Post$alpha_intra[,2] * PlotShade[i])
     # GenericPlot[2,1,i] <- mean(GenericPost)
     # GenericPlot[2,2:3,i] <- HDInterval::hdi(GenericPost)
     # IntraPlot[2,1,i] <- mean(IntraPost)
     # IntraPlot[2,2:3,i] <- HDInterval::hdi(IntraPost)
     # SpecificPlot[2,1,1,i] <- mean(Sp1Post)
     # SpecificPlot[2,1,2:3,i] <- HDInterval::hdi(Sp1Post)
}

AlphaRange <- c(0, 0.02)
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
for(i in 1){
     lines(x = PlotSnow, y = LambdaPlotVals[1,i,,1], col = LambdaCols[i], lwd = 1.5)
     lines(x = PlotSnow, y = LambdaPlotVals[1,i,,2], col = LambdaCols[i], lty = 2)
     lines(x = PlotSnow, y = LambdaPlotVals[1,i,,3], col = LambdaCols[i], lty = 2)
}
#legend("topright", legend = c("DM", "MM"), lty = 1, col = LambdaCols, bty = "n")
text(x = 0.95*SnowRange[1], y = 0.95*LambdaRange[2], 
     labels = expression(bold("a")), cex = 1.5)

# Lambda results for canopy cover
plot(NA, NA, xlim = ShadeRange, ylim = LambdaRange, main = "", xlab = "", ylab = "",
     las = 1, cex.lab = 1.5)
axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
axis(side = 2, at = seq(0, 20, by = 1), labels = FALSE, tcl = -0.25)
for(i in 1:2){
     lines(x = PlotShade, y = LambdaPlotVals[2,i,,1], col = LambdaCols[i], lwd = 1.5)
     lines(x = PlotShade, y = LambdaPlotVals[2,i,,2], col = LambdaCols[i], lty = 2)
     lines(x = PlotShade, y = LambdaPlotVals[2,i,,3], col = LambdaCols[i], lty = 2)
}
legend("topright", legend = c("Bendering", "Perenjori"), lty = 1, col = LambdaCols,
       bty = "n")
text(x = 0.95*ShadeRange[1], y = 0.95*LambdaRange[2], 
     labels = expression(bold("b")), cex = 1.5)

# Alpha results for phosphorous
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
for(i in 1:2){
     lines(x = PlotSnow, y = SpecificPlot[1,i,1,], lwd = 1.5, col = AlphaCols[i])
     #lines(x = PlotSnow, y = SpecificPlot[1,i,2,], lty = 2, col = AlphaCols[i])
     #lines(x = PlotSnow, y = SpecificPlot[1,i,3,], lty = 2, col = AlphaCols[i])
}
legend("topright", legend = SnowLegend, lty = 1, lwd = 1.5, col = c("black", IntraCol, AlphaCols),
       bty = "n")
text(x = 0.95*SnowRange[1], y = 0.95*AlphaRange[2], 
     labels = expression(bold("c")), cex = 1.5)

# Alpha results for canopy cover
plot(NA, NA, xlim = ShadeRange, ylim = AlphaRange, main = "", xlab = "Standardized canopy cover", 
     ylab = "", las = 1, xpd = NA, cex.lab = 1.5)
axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
axis(side = 2, at = seq(0, 0.06, by = 0.005), labels = FALSE, tcl = -0.25)
# Add the generic polygon and line
xCoords <- c(PlotShade, PlotShade[EnvLength:1])
yCoords <- c(GenericPlot[2,2,], GenericPlot[2,3,EnvLength:1])
polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
lines(x = PlotShade, y = GenericPlot[2,1,], lty = 1, lwd = 1.5)
# Add the lines for the intraspecific term
lines(x = PlotShade, y = IntraPlot[2,1,], lwd = 1.5, col = IntraCol)
#lines(x = PlotShade, y = IntraPlot[2,2,], lty = 2, col = IntraCol)
#lines(x = PlotShade, y = IntraPlot[2,3,], lty = 2, col = IntraCol)
# Now the deviating alpha
lines(x = PlotShade, y = SpecificPlot[2,1,1,], lwd = 1.5, col = AlphaCols[1])
#lines(x = PlotShade, y = SpecificPlot[2,1,2,], lty = 2, col = AlphaCols[1])
#lines(x = PlotShade, y = SpecificPlot[2,1,3,], lty = 2, col = AlphaCols[1])
legend("topright", legend = ShadeLegend, lty = 1, lwd = 1.5, col = c("black", IntraCol, AlphaCols),
       bty = "n", inset = -0.02)
text(x = 0.95*ShadeRange[1], y = 0.95*AlphaRange[2], 
     labels = expression(bold("d")), cex = 1.5)
dev.off()

############# Now make a supplemental figure with confidence intervals

SuppFigName <- here("Empirical/WaitziaSuppFig.pdf")
#AlphaRange <- range(c(range(GenericPlot, na.rm = TRUE), range(IntraPlot, na.rm = TRUE), range(SpecificPlot, na.rm = TRUE)))
AlphaRange <- c(0, 0.1)
PhosLetters <- c(expression(bold("a")), expression(bold("c")), expression(bold("e")))
ShadeLetters <- c(expression(bold("b")), expression(bold("d")))
pdf(file = SuppFigName, width = 8, height = 6, onefile = FALSE, paper = "special")
par(mfcol = c(3,2), oma = c(4,4,2,2), mar = c(2,2,2,2))
# Alpha results for phosphorous, divided among three panels for CI visibility
xLabels <- c("", "", "Standardized phosphorous")
for(i in 1:3){
     if(i == 3){
          plot(NA, NA, xlim = PhosRange, ylim = c(0,0.04), main = "", xlab = xLabels[i], 
               ylab = AlphaLab, las = 1, xpd = NA, cex.lab = 1.5)
     }else{
          plot(NA, NA, xlim = PhosRange, ylim = AlphaRange, main = "", xlab = xLabels[i], 
               ylab = AlphaLab, las = 1, xpd = NA, cex.lab = 1.5)
     }
     
     axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
     axis(side = 2, at = seq(0, 0.1, by = 0.005), labels = FALSE, tcl = -0.25)
     # Add the generic polygon and line
     xCoords <- c(PlotPhos, PlotPhos[EnvLength:1])
     yCoords <- c(GenericPlot[1,2,], GenericPlot[1,3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotPhos, y = GenericPlot[1,1,], lty = 1, lwd = 1.5)
     if(i == 1){
          # Add the lines for the intraspecific term
          lines(x = PlotPhos, y = IntraPlot[1,1,], lwd = 1.5, col = IntraCol)
          lines(x = PlotPhos, y = IntraPlot[1,2,], lty = 2, col = IntraCol)
          lines(x = PlotPhos, y = IntraPlot[1,3,], lty = 2, col = IntraCol)
          text(x = 0.5, y = 0.9*AlphaRange[2], labels = PhosLegend[2], col = IntraCol, cex = 1.5)
     }else{
          lines(x = PlotPhos, y = SpecificPlot[1,i-1,1,], lwd = 1.5, col = AlphaCols[i-1])
          lines(x = PlotPhos, y = SpecificPlot[1,i-1,2,], lty = 2, col = AlphaCols[i-1])
          lines(x = PlotPhos, y = SpecificPlot[1,i-1,3,], lty = 2, col = AlphaCols[i-1])
          if(i == 3){
               text(x = 0.5, y = 0.035, labels = PhosLegend[i+1], col = AlphaCols[i-1], cex = 1.5)
               text(x = 0.95*PhosRange[1], y = 0.036, 
                    labels = PhosLetters[i], cex = 1.5)
          }else{
               text(x = 0.5, y = 0.9*AlphaRange[2], labels = PhosLegend[i+1], col = AlphaCols[i-1], cex = 1.5)
          }
     }
     text(x = 0.95*PhosRange[1], y = 0.95*AlphaRange[2], 
          labels = PhosLetters[i], cex = 1.5)
}

# Alpha results for canopy cover, divided among 2 panels
xLabels <- c("", "Standardized canopy cover")
for(i in 1:2){
     plot(NA, NA, xlim = ShadeRange, ylim = AlphaRange, main = "", xlab = xLabels[i], 
          ylab = "", las = 1, xpd = NA, cex.lab = 1.5)
     axis(side = 1, at = seq(-2, 3, by = 0.25), labels = FALSE, tcl = -0.25)
     axis(side = 2, at = seq(0, 0.1, by = 0.005), labels = FALSE, tcl = -0.25)
     # Add the generic polygon and line
     xCoords <- c(PlotShade, PlotShade[EnvLength:1])
     yCoords <- c(GenericPlot[2,2,], GenericPlot[2,3,EnvLength:1])
     polygon(x = xCoords, y = yCoords, border = NA, col = "grey")
     lines(x = PlotShade, y = GenericPlot[2,1,], lty = 1, lwd = 1.5)
     if(i == 1){
          # Add the lines for the intraspecific term
          lines(x = PlotShade, y = IntraPlot[2,1,], lwd = 1.5, col = IntraCol)
          lines(x = PlotShade, y = IntraPlot[2,2,], lty = 2, col = IntraCol)
          lines(x = PlotShade, y = IntraPlot[2,3,], lty = 2, col = IntraCol)
          text(x = 0, y = 0.9*AlphaRange[2], labels = ShadeLegend[2], col = IntraCol, cex = 1.5)
     }else{
          lines(x = PlotShade, y = SpecificPlot[2,1,1,], lwd = 1.5, col = AlphaCols[1])
          lines(x = PlotShade, y = SpecificPlot[2,1,2,], lty = 2, col = AlphaCols[1])
          lines(x = PlotShade, y = SpecificPlot[2,1,3,], lty = 2, col = AlphaCols[1])
          text(x = 0, y = 0.9*AlphaRange[2], labels = ShadeLegend[3], col = AlphaCols[i-1], cex = 1.5)
     }
     text(x = 0.95*ShadeRange[1], y = 0.95*AlphaRange[2], 
          labels = ShadeLetters[i], cex = 1.5)
}
dev.off()


