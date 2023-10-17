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

# Lion demographic dataset
lions.data = read.csv("Data/01_LionsDemographicData.csv")

# Individual capture histories
lions.ch = read.csv("Data/02_LionsCaptureHistories.csv", row.names = 1)
lions.ch = as.matrix(lions.ch)


## 1.4. Covariates ----
# ----------------

pop.size = read.csv("04_Covariate_PopulationSize.csv")$x
rainfall = read.csv("05_Covariate_SeasonalPrecipitation.csv")
rainfall = rainfall$cumul.rain




###########################################################################
#
# 2. Formatting capture histories ----
#
###########################################################################

## 2.1. Getting the first and last capture occasion of each lion ----
# --------------------------------------------------------------

get.first = function(x) min(which(x != 13)) # First observed state
get.last  = function(x){
  
  if(any(x == 11)){which(x == 11)} # Dead recovery occasion
  else{length(x)} # Otherwise, last capture occasion
  
}

lions.first = apply(lions.ch, 1, get.first)
lions.last  = apply(lions.ch, 1, get.last)


## 2.2. Remove lions seen only in the last occasion ----
# -------------------------------------------------

lions.ch = lions.ch[- which(lions.first == lions.last), ]

age = age[- which(lions.first == lions.last), ]
age.unscaled = age.unscaled[- which(lions.first == lions.last), ]
habitat = habitat[- which(lions.first == lions.last), ]

lions.first = apply(lions.ch, 1, get.first)
lions.last  = apply(lions.ch, 1, get.last)


## 2.3. Create group matrix ----
# -------------------------

lions.groups = matrix(NA, nrow = nrow(lions.ch), ncol = ncol(lions.ch), 
                      dimnames = list(rownames(lions.ch), seq(1:60))) # Empty matrix

# Fill in the matrix
for(lion in rownames(lions.groups)){
  
  lions.groups[lion, lions.data$n_census[lions.data$id == lion & lions.data$n_census >= lions_first[lion]]] = 
    lions.data$group[lions.data$id == lion & lions.data$n_census >= lions.first[lion]]
}

# Filling in NAs at t assuming lions stayed in the same group as at t-1
for(lion in rownames(lions.groups)){
  
  for(t in as.numeric(which(is.na(lions.groups[lion, ])))){
    
    if(all(is.na(lions.groups[lion, ]))){break} # If there are only NAs, 
    # go to the next individual
    
    if(t < as.numeric(lions.first[lion])){next} # If the NA is before the 
    # first sighting, continue to 
    # the next iteration
    
    lions.groups[lion, t] = lions.groups[lion, t-1]
    
  }
}

# Remove females only seen as nomad because we cannot assign them group covariates
nomadic.females = which(apply(lions.groups, 1, function(x) all(is.na(x)))) 

lions.ch = lions.ch[- nomadic.females, ]
age = age[- nomadic.females, ]
age.unscaled = age.unscaled[- nomadic.females, ]
habitat = habitat[- nomadic.females, ]
lions_groups = lions_groups[- nomadic.females, ]

# Check that there are no NAs left
lions.groups[which(apply(lions.groups, 1, function(x) any(is.na(x)))), ]
lions.ch[which(apply(lions.groups, 1, function(x) any(is.na(x)))), ]

# Get new first and last sightings
lions.first = apply(lions.ch, 1, get.first)
lions.last  = apply(lions.ch, 1, get.last)




###########################################################################
#
# 3. Loading model output ----
#
###########################################################################

# Get list of output files
files.list = file.info(list.files(path = "", pattern = "LionsFullModel_Output"))[with(file.info(list.files(path = "", pattern = "LionsFullModel_Output")), order(as.POSIXct(mtime))), ]
files.list = rownames(files.list)
files.list = lapply(files.list, FUN = function(x) paste("", x, sep = ""))

# Load files, save the output into an R object, and
# remove the loaded data from the environment to avoid memory issues
load(files.list[[1]])
output1 = res.iter.set
rm(res.iter.set)
load(files.list[[2]])
output2 = res.iter.set
rm(res.iter.set)
load(files.list[[3]])
output3 = res.iter.set
rm(res.iter.set)
load(files.list[[4]])
output4 = res.iter.set
rm(res.iter.set)
load(files.list[[5]])
output5 = res.iter.set
rm(res.iter.set)
load(files.list[[6]])
output6 = res.iter.set
rm(res.iter.set)
load(files.list[[7]])
output7 = res.iter.set
rm(res.iter.set)
load(files.list[[8]])
output8 = res.iter.set
rm(res.iter.set)
load(files.list[[9]])
output9 = res.iter.set
rm(res.iter.set)
load(files.list[[10]])
output10 = res.iter.set
rm(res.iter.set)
load(files.list[[11]])
output11 = res.iter.set
rm(res.iter.set)


# Dividing outputs by chain
output1_chain1 = output1[[1]]
output1_chain2 = output1[[2]]
output1_chain3 = output1[[3]]
output1_chain4 = output1[[4]]

output2_chain1 = output2[[1]]
output2_chain2 = output2[[2]]
output2_chain3 = output2[[3]]
output2_chain4 = output2[[4]]

output3_chain1 = output3[[1]]
output3_chain2 = output3[[2]]
output3_chain3 = output3[[3]]
output3_chain4 = output3[[4]]

output4_chain1 = output4[[1]]
output4_chain2 = output4[[2]]
output4_chain3 = output4[[3]]
output4_chain4 = output4[[4]]

output5_chain1 = output5[[1]]
output5_chain2 = output5[[2]]
output5_chain3 = output5[[3]]
output5_chain4 = output5[[4]]

output6_chain1 = output6[[1]]
output6_chain2 = output6[[2]]
output6_chain3 = output6[[3]]
output6_chain4 = output6[[4]]

output7_chain1 = output7[[1]]
output7_chain2 = output7[[2]]
output7_chain3 = output7[[3]]
output7_chain4 = output7[[4]]

output8_chain1 = output8[[1]]
output8_chain2 = output8[[2]]
output8_chain3 = output8[[3]]
output8_chain4 = output8[[4]]

output9_chain1 = output9[[1]]
output9_chain2 = output9[[2]]
output9_chain3 = output9[[3]]
output9_chain4 = output9[[4]]

output10_chain1 = output10[[1]]
output10_chain2 = output10[[2]]
output10_chain3 = output10[[3]]
output10_chain4 = output10[[4]]

output11_chain1 = output11[[1]]
output11_chain2 = output11[[2]]
output11_chain3 = output11[[3]]
output11_chain4 = output11[[4]]


# Merging datasets
output_chain1 = rbind(output1_chain1, 
                      output2_chain1, 
                      output3_chain1,
                      output4_chain1,
                      output5_chain1,
                      output6_chain1,
                      output7_chain1,
                      output8_chain1,
                      output9_chain1,
                      output10_chain1,
                      output11_chain1)
output_chain2 = rbind(output1_chain2, 
                      output2_chain2, 
                      output3_chain2,
                      output4_chain2,
                      output5_chain2,
                      output6_chain2,
                      output7_chain2,
                      output8_chain2,
                      output9_chain2,
                      output10_chain2,
                      output11_chain2)
output_chain3 = rbind(output1_chain3, 
                      output2_chain3, 
                      output3_chain3,
                      output4_chain3,
                      output5_chain3,
                      output6_chain3,
                      output7_chain3,
                      output8_chain3,
                      output9_chain3,
                      output10_chain3,
                      output11_chain3)
output_chain4 = rbind(output1_chain4, 
                      output2_chain4, 
                      output3_chain4,
                      output4_chain4,
                      output5_chain4,
                      output6_chain4,
                      output7_chain4,
                      output8_chain4,
                      output9_chain4,
                      output10_chain4,
                      output11_chain4)

# Creating an mcmc object with all chains
data.full.mcmc.list = as.mcmc.list(list(as.mcmc(output_chain1),
                                        as.mcmc(output_chain2),
                                        as.mcmc(output_chain3),
                                        as.mcmc(output_chain4)))

# Removing the burn-in phase
data.no.burnin.mcmc.list = as.mcmc.list(list(as.mcmc(output_chain1[10001:nrow(output_chain1),]),
                                             as.mcmc(output_chain2[10001:nrow(output_chain2),]),
                                             as.mcmc(output_chain3[10001:nrow(output_chain3),]),
                                             as.mcmc(output_chain4[10001:nrow(output_chain4),])))

# Full output without burnin as a matrix
lions_output_GLMM = as.matrix(rbind(data.no.burnin.mcmc.list[[1]],
                                    data.no.burnin.mcmc.list[[2]],
                                    data.no.burnin.mcmc.list[[3]],
                                    data.no.burnin.mcmc.list[[4]])) 




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