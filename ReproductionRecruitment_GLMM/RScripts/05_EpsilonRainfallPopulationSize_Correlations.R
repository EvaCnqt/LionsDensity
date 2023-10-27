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


## 1.3. Loading data ----
# ------------------

# MCMC samples
lions_output_fullGLMM = read.csv("Output/ReproductionRecruitmentGLMM_Samples.csv")




###########################################################################
#
# 2. Formatting dataset ----
#
###########################################################################

# Covariates
pop.size = read.csv("Data/02_Covariate_PopulationSize.csv")$x
rainfall = read.csv("Data/03_Covariate_SeasonalPrecipitation.csv")
rainfall = rainfall$cumul.rain




###########################################################################
#
# 3. Plotting settings ----
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
# 4. Calculate correlations between epsilons and rainfall ----
#
###########################################################################

## 4.1. Subset epsilons ----
# ---------------------

epsilons = lions_output_fullGLMM[, grep("epsilon", 
                                        colnames(lions_output_fullGLMM))]

epsilons.repro = epsilons[, grep("epsilon.repro", colnames(epsilons))]
epsilons.rec = epsilons[, grep("epsilon.rec", colnames(epsilons))]


## 4.2. Reproduction probability  ----
# ------------------------------

corr.epsilon.rainfall.pop.size = data.frame(param = "Reproduction probability",
                                            covar = "Rainfall",
                                            corr = apply(epsilons.repro, 1, 
                                                         FUN = function(x) cor(x, rainfall)))


## 4.3. Recruitment  ----
# -----------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Recruitment to\n1 year old",
                                                  covar = "Rainfall",
                                                  corr = apply(epsilons.rec, 1, 
                                                               FUN = function(x) cor(x, rainfall))))




###########################################################################
#
# 5. Calculate correlations between epsilons and population density ----
#
###########################################################################


## 5.1. Reproduction probability  ----
# ------------------------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Reproduction probability",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.repro, 1, 
                                                               FUN = function(x) cor(x, pop.size))))


## 5.2. Recruitment  ----
# -----------------

corr.epsilon.rainfall.pop.size = rbind(corr.epsilon.rainfall.pop.size,
                                       data.frame(param = "Recruitment to\n1 year old",
                                                  covar = "Population size",
                                                  corr = apply(epsilons.rec, 1, 
                                                               FUN = function(x) cor(x, pop.size))))




###########################################################################
#
# 6. Plot correlation coefficients ----
#
###########################################################################

# Factorize parameter names
corr.epsilon.rainfall.pop.size$param = factor(corr.epsilon.rainfall.pop.size$param, 
                                              levels = c("Reproduction probability",
                                                         "Recruitment to\n1 year old"))


# Plot
png(file = "EpsilonPopSizeCorr.png",
    type = "cairo",
    units = "cm",
    width = 8,
    height = 8,
    res = 600,
    bg = "transparent")

ggplot(corr.epsilon.rainfall.pop.size, aes(x = corr)) +
  facet_grid(covar ~ param) +
  geom_density(alpha = 0.2, fill = "#3B0F70FF", col = "#3B0F70FF") +
  geom_vline(xintercept = 0, col = "#FE9F6DFF") +
  xlab("Correlation coefficient") +
  ylab("Density") +
  ylim(0, 7) +
  theme_minimal() %+replace%
  theme(axis.text = element_text(colour = colPlot, size = fontSize, family = font),
        axis.title.x = element_text(colour = colPlot, size = fontSize, 
                                    family = font, 
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = fontSize, 
                                    family = font, 
                                    margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, family = font, 
                                  hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "right", 
        legend.justification = "right", 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(colour = colPlot, size = fontSize, family = font),
        legend.title = element_text(colour = colPlot, size = fontSize, family = font),
        legend.key.size = unit(15, "pt"),
        strip.text = element_text(size = fontSize, margin = margin(b = 2)))

dev.off()
