############################################################################
#
# This script uses samples obtained from chains of an MCMC algorithm.
#
# The aim of this script is to assess correlations between epsilon estimates
# and time-varying covariates (e.g. rainfall, population size).
#
# Author: Eva Conquet
#
###########################################################################

###########################################################################
#
# 1. House keeping and loading libraries and data ----
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading libraries ----
# -----------------------

library(ggplot2)
library(MCMCvis)
library(coda)
library(bayestestR)


## 1.3. Covariates ----
# ----------------

pop.size = read.csv("041_Covariate_PopulationSize.csv")$x
rainfall = read.csv("042_Covariate_SeasonalPrecipitation.csv")
rainfall = rainfall$cumul.rain




###########################################################################
#
# 2. Loading model output ----
#
###########################################################################

lions_output_GLMM = read.csv("Output/MultistateModel_Samples.csv") 




###########################################################################
#
# 4. Plotting settings ----
#
###########################################################################

# Colorblind palette
cbbPalette = c("#000000", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

colBG = "transparent" # Plot background color
colPlot = "black"     # Plot color
font = "Helvetica"    # Plot font
fontSize = 10         # Plot font size

# Define ggplot theme:
theme_general = function(){ 
  
  theme_minimal() %+replace%    # Replace elements we want to change
    
    theme(axis.text = element_text(colour = colPlot, size = 7, family = font),
          axis.title.x = element_text(colour = colPlot, size = 8, family = font, 
                                      margin = margin(t = 15, r = 0, b = 0, l = 0)), 
          axis.title.y = element_text(colour = colPlot, size = 8, family = font, 
                                      margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = colPlot),
          axis.ticks = element_line(colour = colPlot),
          plot.background = element_rect(fill = colBG, color = colBG, size = 0),
          plot.title = element_text(colour = colPlot, size = 10, family = font, 
                                    hjust = 0, 
                                    margin = margin(t = 0, r = 0, b = 15, l = 0)),
          legend.position = c(0.9, 0.9), 
          legend.justification = "right", 
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.text = element_text(colour = colPlot, size = fontSize, family = font),
          legend.title = element_text(colour = colPlot, size = fontSize, family = font),
          legend.key.width = unit(35, "pt"),
          
    ) 
}




###########################################################################
#
# 5. Subset epsilons for each vital rate ----
#
###########################################################################

epsilons = lions_output[, grep("epsilon", colnames(lions_output))]

colnames(epsilons)

epsilons.s.sa1 = epsilons[, grep("epsilon.s.sa1", colnames(epsilons))]
epsilons.s.sa2 = epsilons[, grep("epsilon.s.sa2", colnames(epsilons))]
epsilons.s.af = epsilons[, grep("epsilon.s.af", colnames(epsilons))]
epsilons.s.ym = epsilons[, grep("epsilon.s.ym", colnames(epsilons))]
epsilons.s.nm = epsilons[, grep("epsilon.s.nm", colnames(epsilons))]
epsilons.s.rm = epsilons[, grep("epsilon.s.rm", colnames(epsilons))]
epsilons.s.sa1 = epsilons[, grep("epsilon.s.sa1", colnames(epsilons))]
epsilons.emig.ym = epsilons[, grep("epsilon.emig.ym", colnames(epsilons))]
epsilons.t.ym.nm = epsilons[, grep("epsilon.t.ym.nm", colnames(epsilons))]
epsilons.takeover = epsilons[, grep("epsilon.takeover", colnames(epsilons))]
epsilons.eviction = epsilons[, grep("epsilon.eviction", colnames(epsilons))]
epsilons.dp.pride = epsilons[, grep("epsilon.dp.pride", colnames(epsilons))]
epsilons.dp.nm = epsilons[, grep("epsilon.dp.nm", colnames(epsilons))]
epsilons.dp.dead = epsilons[, grep("epsilon.dp.dead", colnames(epsilons))]




###########################################################################
#
# 6. Calculate correlations between epsilons and rainfall ----
#
###########################################################################

## 6.1. Young subadult survival  ----
# -----------------------------

corr.epsilon.rainfall.pop.size = data.frame(param = "Young-subadult\nsurvival",
                                            covar = "Rainfall",
                                            corr = apply(epsilons.s.sa1, 
                                                         1, 
                                                         FUN = function(x) cor(x, rainfall)))


## 6.2. Old subadult survival  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Old-subadult\nsurvival",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.s.sa2, 
                                                               1,
                                                               FUN = function(x) cor(x, rainfall))))


## 6.3. Adult female survival  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Adult-female\nsurvival",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.s.af, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.4. Young male survival  ----
# -------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young-male\nsurvival",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.s.ym, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.5. Nomadic male survival  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Nomadic-male\nsurvival",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.s.nm,
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.6. Resident male survival  ----
# ----------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Resident-male\nsurvival",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.s.rm, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.7. Young male emigration  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young-male\nemigration",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.emig.ym, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.8. Young male to nomad transition  ----
# ------------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young male to\nnomad\ntransition",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.t.ym.nm, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.9. Nomadic male takeover  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Nomadic-male\ntakeover",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.takeover,
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.10. Resident male eviction  ----
# -----------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Resident-male\neviction",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.eviction, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.11. Pride detection probability  ----
# ----------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Pride detection\nprobability",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.dp.pride,
                                                               1,
                                                               FUN = function(x) cor(x, rainfall))))


## 6.12. Nomadic male detection probability  ----
# -----------------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Nomad\ndetection\nprobability",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.dp.nm, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))


## 6.13. Dead detection probability  ----
# ---------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Dead recovery ",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.dp.dead, 
                                                               1, 
                                                               FUN = function(x) cor(x, rainfall))))




###########################################################################
#
# 7. Calculate correlations between epsilons and population density ----
#
###########################################################################

## 7.1. Young subadult survival  ----
# -----------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young-subadult\nsurvival",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.s.sa1, 
                                                               1,
                                                               FUN = function(x) cor(x, pop.size))))


## 7.2. Old subadult survival  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Old-subadult\nsurvival",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.s.sa2,
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.3. Adult female survival  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Adult-female\nsurvival",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.s.af, 
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.4. Young male survival  ----
# -------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young-male\nsurvival",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.s.ym,
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.5. Nomadic male survival  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Nomadic-male\nsurvival",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.s.nm, 
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.6. Resident male survival  ----
# ----------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Resident-male\nsurvival",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.s.rm, 
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.7. Young male emigration  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young-male\nemigration",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.emig.ym, 
                                                               1,
                                                               FUN = function(x) cor(x, pop.size))))


## 7.8. Young male to nomad transition  ----
# ------------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Young male to\nnomad\ntransition",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.t.ym.nm, 
                                                               1,
                                                               FUN = function(x) cor(x, pop.size))))


## 7.9. Nomadic male takeover  ----
# ---------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Nomadic-male\ntakeover",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.takeover,
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.10. Resident male eviction  ----
# -----------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Resident-male\neviction",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.eviction, 
                                                               1,
                                                               FUN = function(x) cor(x, pop.size))))


## 7.11. Pride detection probability  ----
# ----------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Pride detection\nprobability",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.dp.pride, 
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.12. Nomadic male detection probability  ----
# -----------------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Nomad\ndetection\nprobability",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.dp.nm,
                                                               1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 7.13. Dead detection probability  ----
# ---------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Dead recovery ",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.dp.dead,
                                                               1,
                                                               FUN = function(x) cor(x, pop.size))))




###########################################################################
#
# 8. Plot correlation coefficients ----
#
###########################################################################

# Factorize parameter names
corr.epsilon.rainfall.pop.size$param = factor(corr.epsilon.rainfall.pop.size$param, 
                                              levels = unique(corr.epsilon.rainfall.pop.size$param))


# Plot
png(file = "EpsilonPopSizeCorr.png",
    type = "cairo",
    units = "cm",
    width = 23,
    height = 6,
    res = 600,
    bg = "transparent")

ggplot(corr.epsilon.rainfall.pop.size, aes(x = corr)) +
  facet_grid(covar ~ param) +
  geom_density(alpha = 0.2, fill = "#3B0F70FF", col = "#3B0F70FF") +
  geom_vline(xintercept = 0, col = "#FE9F6DFF") +
  xlab("Correlation coefficient") +
  ylab("Density") +
  theme_minimal() %+replace%    
  
  theme(axis.text = element_text(colour = colPlot, size = fontSize, family = font),
        axis.title.x = element_text(colour = colPlot, size = fontSize, family = font, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = fontSize, family = font, margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, family = font, hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "right", 
        legend.justification = "right", 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(colour = colPlot, size = fontSize, family = font),
        legend.title = element_text(colour = colPlot, size = fontSize, family = font),
        legend.key.size = unit(15, "pt"),
        strip.text = element_text(size = fontSize, margin = margin(b = 2)))

dev.off()