############################################################################
#
# This script uses samples obtained from chains of an MCMC algorithm.
#
# The aim of this script is to check the convergence of the chains and obtain 
# information on parameters posterior distributions.
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

rm(list = ls(all.names = T))


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

# Year and season
year = read.csv("Data/031_Covariate_Year.csv", stringsAsFactors = F)$x
year = rep(c(1:30), each = 2)

season = read.csv("Data/032_Covariate_Season.csv", stringsAsFactors = F)$x


# Habitat
habitat = read.csv("Data/033_Covariate_Habitat.csv", row.names = 1)
habitat = as.matrix(habitat)


# Density-dependent covariates

# Number of females in a pride
nb.af.pride.unscaled = read.csv("Data/034_Covariate_NbAFpride.csv", 
                                stringsAsFactors = F, 
                                row.names = 1)
nb.af.pride.unscaled = as.matrix(nb.af.pride.unscaled)
range(nb.af.pride.unscaled, na.rm = T)
nb.af.pride = (nb.af.pride.unscaled - mean(nb.af.pride.unscaled, na.rm = T)) / 
  (2 * sd(nb.af.pride.unscaled, na.rm = T)) # Standardize covariate


# Age
age.unscaled = read.csv("Data/035_Covariate_Age.csv", row.names = 1)
age.unscaled = as.matrix(age.unscaled)
range(age.unscaled, na.rm = T)
age = (age.unscaled - mean(age.unscaled, na.rm = T)) / 
  (2 * sd(age.unscaled, na.rm = T)) # Standardize covariate


# Male coalition size 
coal.size.unscaled = read.csv("Data/036_Covariate_CoalSize.csv", row.names = 1)
coal.size.unscaled = as.matrix(coal.size.unscaled)
range(coal.size.unscaled, na.rm = T)
coal.size = (coal.size.unscaled - mean(coal.size.unscaled, na.rm = T)) / 
  (2 * sd(coal.size.unscaled, na.rm = T)) # Standardize covariate


# Number of nomadic coalitions in the home range of a pride or 
# a resident male coalition
nb.nm.coal.hr.unscaled = read.csv("Data/037_Covariate_NbNMCoalHR.csv", row.names = 1)
nb.nm.coal.hr.unscaled = as.matrix(nb.nm.coal.hr.unscaled)
range(nb.nm.coal.hr.unscaled, na.rm = T)
nb.nm.coal.hr = (nb.nm.coal.hr.unscaled - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / 
  (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)) # Standardize covariate




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
# 4. Checking for convergence by plotting all three chains ----
# for each parameter
#
###########################################################################

# List of parameters
parameters.lionKFDHMM = c(# Intercepts
                          "mu.s.sa1",    # Young-subadult survival
                          "mu.s.sa2f",   # Female old-subadult survival 
                          "mu.s.sa2m",   # Male old-subadult survival
                          "mu.s.af",     # Adult-female survival
                          "mu.s.ym",     # Young-male survival
                          "mu.s.nm",     # Nomadic-male survival
                          "mu.s.rm",     # Resident-male survival
                          "mu.emig.ym",  # Young-male emigration
                          "mu.t.ym.nm",  # Young-male transition to nomadic male
                          "mu.takeover", # Nomadic-male takeover
                          "mu.eviction", # Resident-male eviction 
                          "mu.dp.pride", # Detection of pride individual
                          "mu.dp.nm",    # Detection of nomadic male
                          "mu.dp.dead",  # Detection of dead individual
                          
                          # Betas
                          # Young-subadult survival
                          "s.sa1.beta.nb.nm.coal.hr",       # Number of nomadic coalitions in the home range
                          "s.sa1.beta.nb.af.pride",         # Number of adult females in the pride
                          "s.sa1.beta.habitat.woodland",    # Habitat (NA for grassland which is the reference level)
                          # Old-subadult survival
                          "s.sa2.beta.nb.nm.coal.hr",       # Number of nomadic coalitions in the home range 
                          "s.sa2.beta.nb.af.pride",         # Number of adult females in the pride
                          "s.sa2.beta.habitat.woodland",    # Habitat (NA for grassland which is the reference level)
                          # Adult-female survival
                          "s.af.beta.age",                  # Age
                          "s.af.beta.nb.af.pride",          # Number of adult females in the pride
                          "s.af.beta.nb.nm.coal.hr",        # Number of nomadic coalitions in the home range
                          "s.af.beta.habitat.woodland",     # Habitat (NA for grassland which is the reference level)
                          # Young-male survival
                          "s.ym.beta.nb.nm.coal.hr",        # Number of nomadic coalitions in the home range
                          "s.ym.beta.nb.af.pride",          # Number of adult females in the pride
                          "s.ym.beta.habitat.woodland",     # Habitat (NA for grassland which is the reference level)
                          # Nomadic-male survival
                          "s.nm.beta.coal.size",            # Coalition size
                          "s.nm.beta.habitat.woodland",     # Habitat (NA for grassland which is the reference level)
                          # Resident-male survival
                          "s.rm.beta.coal.size",            # Coalition size
                          "s.rm.beta.nb.nm.coal.hr",        # Number of nomadic coalitions in the home range
                          "s.rm.beta.habitat.woodland",     # Habitat (NA for grassland which is the reference level)
                          # Young-male emigration
                          "emig.ym.beta.habitat.woodland",  # Habitat (NA for grassland which is the reference level)
                          # Nomadic-male takeover
                          "takeover.beta.coal.size",        # Coalition size
                          "takeover.beta.habitat.woodland", # Habitat (NA for grassland which is the reference level) 
                          # Resident-male eviction
                          "eviction.beta.coal.size",        # Coalition size
                          "eviction.beta.nb.nm.coal.hr",    # Number of nomadic coalitions in the home range
                          "eviction.beta.habitat.woodland", # Habitat (NA for grassland which is the reference level) 
                          # Pride-individual detection
                          "dp.pride.beta.habitat.woodland", # Habitat (NA for grassland which is the reference level) 
                          # Nomadic-male detection
                          "dp.nm.beta.habitat.woodland",    # Habitat (NA for grassland which is the reference level)
                          
                          # Random effect sigmas
                          "sigma.s.sa1",    # Young-subadult survival
                          "sigma.s.sa2",    # Old-subadult survival
                          "sigma.s.af",     # Adult-female survival
                          "sigma.s.ym",     # Young-male survival
                          "sigma.s.nm",     # Nomadic-male survival
                          "sigma.s.rm",     # Resident-male survival
                          "sigma.emig.ym",  # Young-male emigration
                          "sigma.t.ym.nm",  # Young-male transition to nomadic male
                          "sigma.takeover", # Nomadic-male takeover
                          "sigma.eviction", # Resident-male eviction
                          "sigma.dp.pride", # Pride-individual detection
                          "sigma.dp.nm",    # Nomadic-male detection
                          "sigma.dp.dead",  # Dead-individual detection
                          
                          # Random effect epsilons
                          "epsilon.s.sa1",    # Young-subadult survival  
                          "epsilon.s.sa2",    # Old-subadult survival
                          "epsilon.s.af",     # Adult-female survival
                          "epsilon.s.ym",     # Young-male survival
                          "epsilon.s.nm",     # Nomadic-male survival
                          "epsilon.s.rm",     # Resident-male survival
                          "epsilon.emig.ym",  # Young-male emigration
                          "epsilon.t.ym.nm",  # Young-male transition to nomadic male
                          "epsilon.takeover", # Nomadic-male takeover
                          "epsilon.eviction", # Resident-male eviction
                          "epsilon.dp.pride", # Pride-individual detection
                          "epsilon.dp.nm",    # Nomadic-male detection
                          "epsilon.dp.dead")  # Dead-individual detection


# Remove epsilons and sigmas for the traceplots
params.to.plot = parameters.lionKFDHMM[- grep("epsilon", parameters.lionKFDHMM)]
params.to.plot = params.to.plot[- grep("sigma", params.to.plot)]


# Traceplots of all chains
MCMCtrace(data.full.mcmc.list, 
          iter = nrow(output_chain1), # Plot all iterations 
          params = params.to.plot, # Parameters to plot
          priors = runif(nrow(output_chain1) * 4, -10, 10), # Priors to get the prior-posterior overlaps 
          filename = "CMR_MCMCtrace.pdf")




###########################################################################
#
# 5. Plotting settings ----
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

# Define ggplot theme
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
# 6. Prior-posterior overlap ----
#
###########################################################################

# Empty dataframe to store the prior and posterior data
prior_post_param_df = data.frame(parameter = NA,
                                 season = NA,
                                 prior_post = NA,
                                 distribution = NA)[0, ]

# Parameters name to plot
param_full_name = data.frame(param_short = params.to.plot,
                             param_full = c("Mean young-subadult\nsurvival",
                                            "Mean female\nold-subadult survival",
                                            "Mean male\nold-subadult survival",
                                            "Mean adult-female\nsurvival",
                                            "Mean young-male survival",
                                            "Mean nomadic-male\nsurvival",
                                            "Mean resident-male\nsurvival",
                                            "Mean young-male\nemigration",
                                            "Mean transition from young\nto nomadic male",
                                            "Mean nomadic-male\ntakeover",
                                            "Mean resident-male\neviction",
                                            "Mean pride\ndetection probability",
                                            "Mean nomadic-male\ndetection probability",
                                            "Mean dead\nrecovery probability",
                                            "Young-subadult survival:\nNb nomadic coalitions HR",
                                            "Young-subadult survival:\nNb females pride",
                                            "Young-subadult survival:\nHabitat (woodland)",
                                            "Old-subadult survival:\nNb nomadic coalitions HR",
                                            "Old-subadult survival:\nNb females pride",
                                            "Old-subadult survival:\nHabitat (woodland)",
                                            "Adult-female survival:\nAge",
                                            "Adult-female survival:\nNb females pride",
                                            "Adult-female survival:\nNb nomadic coalitions HR",
                                            "Adult-female survival:\nHabitat (woodland)",
                                            "Young-male survival:\nNb nomadic coalitions HR",
                                            "Young-male survival:\nNb females pride",
                                            "Young-male survival:\nHabitat (woodland)",
                                            "Nomadic-male survival:\nCoalition size",
                                            "Nomadic-male survival:\nHabitat (woodland)",
                                            "Resident-male survival:\nCoalition size",
                                            "Resident-male survival:\nNb nomadic coalitions HR",
                                            "Resident-male survival:\nHabitat (woodland)",
                                            "Young-male emigration:\nHabitat (woodland)",
                                            "Nomadic-male takeover:\nCoalition size",
                                            "Nomadic-male takeover:\nHabitat (woodland)",
                                            "Resident-male eviction:\nCoalition size",
                                            "Resident-male eviction:\nNb nomadic coalitions HR",
                                            "Resident-male eviction:\nHabitat (woodland)",
                                            "Pride detection:\nHabitat (woodland)",
                                            "Nomadic-male detection:\nHabitat (woodland)"))


# Fill in dataframe with prior and posterior distributions
for(param in params.to.plot){
  
  # Full parameter name
  full_param = param_full_name$param_full[which(param_full_name$param_short == param)]
  
  # For the parameter of the effect of nomadic males on old subadults,
  # which has only one average estimate
  if(param == "s.sa2.beta.nb.nm.coal.hr"){
    
    prior_post_param = data.frame(parameter = full_param,
                                  season = NA,
                                  prior_post = rep(c("prior", "posterior"), 
                                                   each = nrow(lions_output_GLMM)),
                                  distribution = c(runif(nrow(lions_output_GLMM), -10, 10), 
                                                   lions_output_GLMM[, grep(param, colnames(lions_output_GLMM))]))
    
  }
  
  # For all other parameters, which have one estimate per season
  else{
    
    prior_post_param = data.frame(parameter = full_param,
                                  season = rep(c("wet", "dry"), 
                                               each = 2 * nrow(lions_output_GLMM)),
                                  prior_post = rep(rep(c("prior", "posterior"), 
                                                       each = nrow(lions_output_GLMM)), 2),
                                  distribution = c(runif(nrow(lions_output_GLMM), -10, 10), 
                                                   lions_output_GLMM[, grep(param, colnames(lions_output_GLMM))][, 1], 
                                                   runif(nrow(lions_output_GLMM), -10, 10), 
                                                   lions_output_GLMM[, grep(param, colnames(lions_output_GLMM))][, 2]))

  }
  
  # Merge data to the full dataframe
  prior_post_param_df = rbind(prior_post_param_df, prior_post_param)
}


# Format prior-posterior dataframe
prior_post_param_df$season = paste0("\n(", prior_post_param_df$season, " season)")
prior_post_param_df$param_seas = paste(prior_post_param_df$parameter, prior_post_param_df$season, sep = " ")
prior_post_param_df$param_seas = factor(prior_post_param_df$param_seas, levels = unique(prior_post_param_df$param_seas))


# Prior-posterior overlap plots
# (the percentages have been aded manually using the MCMCtrace plot information)

temp = seq(1, 40, 6) # Limits of numbers of plots per page

for(i in 1:length(temp)){
  
  j = ifelse(i != length(temp), temp[i] + 5, 
             length(unique(prior_post_param_df$parameter))) # Index of last plot on the current page

  # Prior-posterior overlap plot
  prior_post_plot = 
  ggplot(prior_post_param_df[which(prior_post_param_df$parameter %in% 
                                   unique(prior_post_param_df$parameter)[(temp[i]):(j)]), ], 
                             aes(x = distribution, 
                                 fill = prior_post, colour = prior_post)) +
  facet_wrap(~ param_seas, ncol = 3, scale = "free_y") +
  geom_density(alpha = 0.2) +
  scale_fill_manual(name = "Distribution",
                    labels = c("Posterior", "Prior"),
                    values = c("#FE9F6DFF", "#3B0F70FF")) +
  scale_colour_manual(name = "Distribution",
                      labels = c("Posterior", "Prior"),
                      values = c("#FE9F6DFF", "#3B0F70FF")) +
  scale_y_continuous(n.breaks = 4) +
  xlab("Parameter value") +
  ylab("Density") +
  theme_minimal() %+replace%    
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
        plot.title = element_text(colour = colPlot, size = 8, family = font, 
                                  hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "right", 
        legend.justification = "right",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(colour = colPlot, size = 7, family = font),
        legend.title = element_text(colour = colPlot, size = 8, family = font),
        legend.key.size = unit(15, "pt"),
        strip.text = element_text(size = 8, margin = margin(b = 10)))


png(filename = paste0("PriorPosteriorOverlap_", i, ".png"), 
    width = 16, 
    height = 12, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

print(prior_post_plot)

dev.off()

}



###########################################################################
#
# 7. Effect size and prediction plots ----
#
###########################################################################

## 7.1. Young subadult survival ----
# -----------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_s.sa1 = lions_output_GLMM[, grep("s.sa1.", colnames(lions_output_GLMM))]
lions_output_s.sa1.epsilons = lions_output_s.sa1[, grep("epsilon", colnames(lions_output_s.sa1))]
lions_output_s.sa1 = lions_output_s.sa1[, - grep("epsilon", colnames(lions_output_s.sa1))]
lions_output_s.sa1 = lions_output_s.sa1[, - grep("sigma", colnames(lions_output_s.sa1))]


# Mean, median,  and 90% credible interval of the estimates
lions_output_s.sa1_estimates = data.frame(parameter = colnames(lions_output_s.sa1), 
                                          mean = apply(lions_output_s.sa1, 2, mean),
                                          median = apply(lions_output_s.sa1, 2, median),
                                          lowerCI = apply(lions_output_s.sa1, 2, 
                                                          FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                          upperCI = apply(lions_output_s.sa1, 2, 
                                                          FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_s.sa1_estimates$parameter_clean = c("Mean survival\n(wet season in the grassland)",
                                                 "Mean survival\n(dry season in the grassland)",
                                                 "Habitat\n(grassland, wet season)",
                                                 "Habitat\n(grassland, dry season)",
                                                 "Habitat\n(woodland, wet season)",
                                                 "Habitat\n(woodland, dry season)",
                                                 "Number of females\nin the pride (wet season)",
                                                 "Number of females\nin the pride (dry season)",
                                                 "Number of nomadic\ncoalitions in the home range (wet season)",
                                                 "Number of nomadic\ncoalitions in the home range (dry season)")


# Add season and parameter category
lions_output_s.sa1_estimates$season = rep(c("Wet", "Dry"), 5)
lions_output_s.sa1_estimates$parameter_plot = rep(c("Mean survival", 
                                                    "Habitat (grassland)", 
                                                    "Habitat\n(woodland)", 
                                                    "Number of females\nin the pride", 
                                                    "Number of nomadic\ncoalitions in the home range"), 
                                                  each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_s.sa1_estimates$parameter_plot))),
                               variable = rep(c("Mean survival", "Habitat (grassland)", 
                                                "Habitat\n(woodland)", "Number of females\nin the pride", 
                                                "Number of nomadic\ncoalitions in the home range"), 
                                              each = 2 * n),
                               posterior = c(lions_output_s.sa1[1:nrow(lions_output_s.sa1), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean survival", 
                                               "Number of females\nin the pride", 
                                               "Number of nomadic\ncoalitions in the home range", 
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.1.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_s.sa1_estimates$parameter_plot,
                     season = lions_output_s.sa1_estimates$season,
                     mean = lions_output_s.sa1_estimates$mean,
                     median = lions_output_s.sa1_estimates$median,
                     BCI90_upper = lions_output_s.sa1_estimates$upperCI,
                     BCI90_lower = lions_output_s.sa1_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean survival", 
                                     "Number of females\nin the pride", 
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "SA1_Survival.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

sa1surv_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, 
                           fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) + 
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, aes(x = median, y = variable, # Median
                                 color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping zero",
                     labels = c("No", "Yes"),
                     values = c(1, 0.4)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 2), limits = c(-5.7, 10)) +
  theme_general() +
  ggtitle("Young-subadult survival") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

sa1surv_plot

dev.off()


## 7.1.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
ys.survival.estimate = function(season = "wet",
                                habitat = "grassland",
                                nb.af.pride = 0,
                                nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.sa1", 
                                         colnames(lions_output_GLMM))][, 1]
    
    beta.nb.af.pride = lions_output_GLMM[, grep("s.sa1.beta.nb.af", 
                                                colnames(lions_output_GLMM))][, 1]
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.sa1.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_s.sa1.epsilons[, seq(1, 59, 2)], 1, median)

    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.sa1.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.sa1", 
                                         colnames(lions_output_GLMM))][, 2]
    
    beta.nb.af.pride = lions_output_GLMM[, grep("s.sa1.beta.nb.af", 
                                                colnames(lions_output_GLMM))][, 2]
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.sa1.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_s.sa1.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.sa1.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred = mean.surv + beta.habitat + beta.nb.af.pride * nb.af.pride + 
         beta.nb.nm.coal.hr * nb.nm.coal.hr + median.epsilons
  
  return(pred)
  
}


# Average survival per habitat and season

# Empty prediction dataframe
ys.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                          habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
ys.season.habitat.estimates$pred = plogis(apply(apply(ys.season.habitat.estimates, 
                                                      1, 
                                                      FUN = function(x) ys.survival.estimate(season = x[1], habitat = x[2])), 
                                                2, median))
ys.season.habitat.estimates$lwr = plogis(apply(apply(ys.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) ys.survival.estimate(season = x[1], habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
ys.season.habitat.estimates$upr = plogis(apply(apply(ys.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) ys.survival.estimate(season = x[1], habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_SA1_Survival_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ys.season.habitat.estimates, aes(x = season, y = pred, 
                                        colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Season") +
  ylab("Young-subadult\nsurvival probability") +
  labs(title = "Young-subadult survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


# Average survival with number of adult females in the pride

# Empty prediction dataframe
ys.nb.af.estimates = expand.grid(season = c("wet", "dry"),
                                 nb.af.pride = (seq(1, 12) - mean(nb.af.pride.unscaled, na.rm = T)) / (2 * sd(nb.af.pride.unscaled, na.rm = T)),
                                 pred = NA)

# Fill in prediction and credible intervals
ys.nb.af.estimates$pred = plogis(apply(apply(ys.nb.af.estimates, 
                                             1, 
                                             FUN = function(x) ys.survival.estimate(season = x[1], nb.af.pride = as.numeric(x[2]))), 
                                       2, median))
ys.nb.af.estimates$lwr = plogis(apply(apply(ys.nb.af.estimates,
                                            1, 
                                            FUN = function(x) ys.survival.estimate(season = x[1], nb.af.pride = as.numeric(x[2]))), 
                                      2, FUN = function(x) quantile(x, probs = 0.05)))
ys.nb.af.estimates$upr = plogis(apply(apply(ys.nb.af.estimates, 
                                            1, 
                                            FUN = function(x) ys.survival.estimate(season = x[1], nb.af.pride = as.numeric(x[2]))), 
                                      2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_SA1_Survival_NbAFPride.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ys.nb.af.estimates, aes(x = nb.af.pride, y = pred, 
                               colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, 
                  colour = season, fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = unique(ys.nb.af.estimates$nb.af.pride),
                     labels = seq(1, 12)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Number of females\nin the pride") +
  ylab("Young-subadult\nsurvival probability") +
  labs(title = "Young-subadult survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none") 

dev.off()  


# Average survival with number of nomadic coalitions in the home range

# Empty prediction dataframe
ys.nb.nm.coal.estimates = expand.grid(season = c("wet", "dry"),
                                      nb.nm.coal.hr = (seq(0, 6) - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)),
                                      pred = NA)

# Fill in prediction and credible intervals
ys.nb.nm.coal.estimates$pred = plogis(apply(apply(ys.nb.nm.coal.estimates,
                                                  1, 
                                                  FUN = function(x) ys.survival.estimate(season = x[1], nb.nm.coal.hr = as.numeric(x[2]))), 
                                            2, median))
ys.nb.nm.coal.estimates$lwr = plogis(apply(apply(ys.nb.nm.coal.estimates, 
                                                 1, 
                                                 FUN = function(x) ys.survival.estimate(season = x[1], nb.nm.coal.hr = as.numeric(x[2]))),
                                           2, FUN = function(x) quantile(x, probs = 0.05)))
ys.nb.nm.coal.estimates$upr = plogis(apply(apply(ys.nb.nm.coal.estimates,
                                                 1, 
                                                 FUN = function(x) ys.survival.estimate(season = x[1], nb.nm.coal.hr = as.numeric(x[2]))), 
                                           2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_SA1_Survival_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ys.nb.nm.coal.estimates, aes(x = nb.nm.coal.hr, y = pred, 
                                    colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, 
                  fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = unique(ys.nb.nm.coal.estimates$nb.nm.coal.hr),
                     labels = seq(0, 6)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Young-subadult\nsurvival probability") +
  labs(title = "Young-subadult survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


## 7.2. Old subadult survival ----
# -----------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_s.sa2 = lions_output_GLMM[, grep("s.sa2.", 
                                              colnames(lions_output_GLMM))]
lions_output_s.sa2.epsilons = lions_output_s.sa2[, grep("epsilon", 
                                                        colnames(lions_output_s.sa2))]
lions_output_s.sa2 = lions_output_s.sa2[, - grep("epsilon", 
                                                 colnames(lions_output_s.sa2))]
lions_output_s.sa2 = lions_output_s.sa2[, - grep("sigma", 
                                                 colnames(lions_output_s.sa2))]

# Mean, median,  and 90% credible interval of the estimates
lions_output_s.sa2_estimates = data.frame(parameter = colnames(lions_output_s.sa2), 
                                          mean = apply(lions_output_s.sa2, 2, mean),
                                          median = apply(lions_output_s.sa2, 2, median),
                                          lowerCI = apply(lions_output_s.sa2, 2, 
                                                          FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                          upperCI = apply(lions_output_s.sa2, 2, 
                                                          FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_s.sa2_estimates$parameter_clean = c("Mean female survival\n(wet season in the grassland)",
                                                 "Mean female survival\n(dry season in the grassland)",
                                                 "Mean male survival\n(wet season in the grassland)",
                                                 "Mean male survival\n(dry season in the grassland)",
                                                 "Habitat\n(grassland, wet season)",
                                                 "Habitat\n(grassland, dry season)",
                                                 "Habitat\n(woodland, wet season)",
                                                 "Habitat\n(woodland, dry season)",
                                                 "Number of females\nin the pride (wet season)",
                                                 "Number of females\nin the pride (dry season)",
                                                 "Number of nomadic\ncoalitions in the home range")

# Add season and parameter category
lions_output_s.sa2_estimates$season = c(rep(c("Wet", "Dry"), 5), NA)
lions_output_s.sa2_estimates$parameter_plot = c(rep(c("Mean female survival", 
                                                      "Mean male survival", 
                                                      "Habitat (grassland)", 
                                                      "Habitat (woodland)", 
                                                      "Number of females\nin the pride"), 
                                                    each = 2), 
                                                "Number of nomadic\ncoalitions in the home range")

# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = c(rep(rep(c("Wet", "Dry"), each = n), 
                                              length(unique(lions_output_s.sa2_estimates$parameter_plot)) - 1), 
                                          rep(NA, n)),
                               variable = c(rep(c("Mean female survival", 
                                                  "Mean male survival",
                                                  "Habitat (grassland)",
                                                  "Habitat (woodland)", 
                                                  "Number of females\nin the pride"),
                                                each = 2 * n), 
                                            rep("Number of nomadic\ncoalitions in the home range", n)),
                               posterior = c(lions_output_s.sa2[1:nrow(lions_output_s.sa2), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean female survival", 
                                               "Mean male survival", 
                                               "Number of females\nin the pride",
                                               "Number of nomadic\ncoalitions in the home range", 
                                               "Habitat (woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.2.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_s.sa2_estimates$parameter_plot,
                     season = lions_output_s.sa2_estimates$season,
                     mean = lions_output_s.sa2_estimates$mean,
                     median = lions_output_s.sa2_estimates$median,
                     BCI90_upper = lions_output_s.sa2_estimates$upperCI,
                     BCI90_lower = lions_output_s.sa2_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean female survival", "Mean male survival", 
                                     "Number of females\nin the pride", 
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Habitat (woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}

# Plot effect size
png(filename = "SA2_Survival.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

sa2surv_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plots
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) + 
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, aes(x = median, y = variable, # Median
                                 color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping zero",
                     labels = c("Yes", "No"),
                     values = c(1, 0.4)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 1), limits = c(-4.6, 10)) +
  theme_general() +
  ggtitle("Old-subadult survival") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

sa2surv_plot

dev.off()


## 7.2.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
os.survival.estimate = function(sex = "female",
                                season = "wet",
                                habitat = "grassland",
                                nb.af.pride = 0,
                                nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get mean female survival
    if(sex == "female"){
      
      mean.surv = lions_output_GLMM[, grep("mu.s.sa2f", 
                                           colnames(lions_output_GLMM))][, 1]
      
    }
    
    # Get mean male survival
    else{
      
      mean.surv = lions_output_GLMM[, grep("mu.s.sa2m", 
                                           colnames(lions_output_GLMM))][, 1]
      
    }
    
    # Get wet-season parameter values
    beta.nb.af.pride = lions_output_GLMM[, grep("s.sa2.beta.nb.af", 
                                                colnames(lions_output_GLMM))][, 1]
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.sa2.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_s.sa2.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.sa2.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get mean female survival
    if(sex == "female"){
      
      mean.surv = lions_output_GLMM[, grep("mu.s.sa2f", 
                                           colnames(lions_output_GLMM))][, 2]
      
    }
    
    # Get mean male survival
    else{
      
      mean.surv = lions_output_GLMM[, grep("mu.s.sa2m", 
                                           colnames(lions_output_GLMM))][, 2]
      
    }
    
    # Get dry-season parameter values
    beta.nb.af.pride = lions_output_GLMM[, grep("s.sa2.beta.nb.af", 
                                                colnames(lions_output_GLMM))][, 2]
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.sa2.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_s.sa2.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dty-season habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.sa2.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred = mean.surv + beta.habitat + beta.nb.af.pride * nb.af.pride + 
         beta.nb.nm.coal.hr * nb.nm.coal.hr + median.epsilons
  
  return(pred)
  
}


# Average survival per habitat and season

# Empty prediction dataframe
os.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                          habitat = c("grassland", "woodland"),
                                          sex = c("female", "male"))

# Fill in prediction and credible intervals
os.season.habitat.estimates$pred = plogis(apply(apply(os.season.habitat.estimates, 
                                                      1, 
                                                      FUN = function(x) os.survival.estimate(season = x[1], 
                                                                                             habitat = x[2],
                                                                                             sex = x[3])), 
                                                2, median))
os.season.habitat.estimates$lwr = plogis(apply(apply(os.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) os.survival.estimate(season = x[1],
                                                                                            habitat = x[2],
                                                                                            sex = x[3])), 
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
os.season.habitat.estimates$upr = plogis(apply(apply(os.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) os.survival.estimate(season = x[1],
                                                                                            habitat = x[2], 
                                                                                            sex = x[3])), 
                                               2, FUN = function(x) quantile(x, probs = 0.95)))

# Aggregate over sexes
os.season.habitat.estimates = aggregate(cbind(pred, lwr, upr) ~ habitat + season, 
                                        data = os.season.habitat.estimates, mean)


# Prediction plot
png(filename = "Predictions_SA2_Survival_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(os.season.habitat.estimates, aes(x = season, y = pred, 
                                        colour = habitat)) +
  facet_wrap(~ sex, labeller = labeller(sex = c("female" = "Female",
                                                "male" = "Male"))) + 
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Old-subadult\nsurvival probability") +
  labs(title = "Old-subadult survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average survival with number of adult females in the pride

# Empty prediction dataframe
os.nb.af.estimates = expand.grid(season = c("wet", "dry"),
                                 nb.af.pride = (seq(1, 12) - mean(nb.af.pride.unscaled, na.rm = T)) / (2 * sd(nb.af.pride.unscaled, na.rm = T)),
                                 sex = c("female", "male"),
                                 pred = NA)

# Fill in prediction and credible intervals
os.nb.af.estimates$pred = plogis(apply(apply(os.nb.af.estimates, 
                                             1, 
                                             FUN = function(x) os.survival.estimate(season = x[1], 
                                                                                    nb.af.pride = as.numeric(x[2]), 
                                                                                    sex = x[3])), 
                                       2, median))
os.nb.af.estimates$lwr = plogis(apply(apply(os.nb.af.estimates, 
                                            1, 
                                            FUN = function(x) os.survival.estimate(season = x[1],
                                                                                   nb.af.pride = as.numeric(x[2]), 
                                                                                   sex = x[3])), 
                                      2, FUN = function(x) quantile(x, probs = 0.05)))
os.nb.af.estimates$upr = plogis(apply(apply(os.nb.af.estimates, 
                                            1, 
                                            FUN = function(x) os.survival.estimate(season = x[1], 
                                                                                   nb.af.pride = as.numeric(x[2]), 
                                                                                   sex = x[3])), 
                                      2, FUN = function(x) quantile(x, probs = 0.95)))

# Aggregate over sexes
os.nb.af.estimates = aggregate(cbind(pred, lwr, upr) ~ season + nb.af.pride, 
                               data = os.nb.af.estimates, FUN = mean)

# Prediction plot
png(filename = "Predictions_SA2_Survival_NbAFPride.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(os.nb.af.estimates, aes(x = nb.af.pride, y = pred, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(os.nb.af.estimates$nb.af.pride),
                     labels = seq(1, 12)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Number of females\nin the pride") +
  ylab("Old-subadult\nsurvival probability") +
  labs(title = "Old-subadult survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


# Average survival with number of nomadic coalitions in the home range

# Empty prediction dataframe
os.nb.nm.coal.estimates = expand.grid(season = c("wet", "dry"),
                                      nb.nm.coal.hr = (seq(0, 6) - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)),
                                      sex = c("female", "male"),
                                      pred = NA)

# Fill in prediction and credible intervals
os.nb.nm.coal.estimates$pred = plogis(apply(apply(os.nb.nm.coal.estimates,
                                                  1, 
                                                  FUN = function(x) os.survival.estimate(season = x[1],
                                                                                         nb.nm.coal.hr = as.numeric(x[2]), 
                                                                                         sex = x[3])), 
                                            2, median))
os.nb.nm.coal.estimates$lwr = plogis(apply(apply(os.nb.nm.coal.estimates, 
                                                 1,
                                                 FUN = function(x) os.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]),
                                                                                        sex = x[3])), 
                                           2, FUN = function(x) quantile(x, probs = 0.05)))
os.nb.nm.coal.estimates$upr = plogis(apply(apply(os.nb.nm.coal.estimates, 
                                                 1, 
                                                 FUN = function(x) os.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]),
                                                                                        sex = x[3])), 
                                           2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_SA2_Survival_NbNMCoalHR.png", 
    width = 3000, 
    height = 1500, 
    units = "px", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(os.nb.nm.coal.estimates, aes(x = nb.nm.coal.hr, y = pred, colour = season)) +
  facet_wrap(~ sex, labeller = labeller(sex = c("female" = "Female",
                                                "male" = "Male"))) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(os.nb.nm.coal.estimates$nb.nm.coal.hr),
                     labels = seq(0, 6)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Old-subadult\nsurvival probability") +
  labs(title = "Old-subadult survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


## 7.3. Adult female survival ----
# ---------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_s.af = lions_output_GLMM[, grep("s.af.", 
                                             colnames(lions_output_GLMM))]
lions_output_s.af.epsilons = lions_output_s.af[, grep("epsilon", 
                                                      colnames(lions_output_s.af))]

lions_output_s.af = lions_output_s.af[, - grep("epsilon", 
                                               colnames(lions_output_s.af))]
lions_output_s.af = lions_output_s.af[, - grep("sigma", 
                                               colnames(lions_output_s.af))]

# Mean, median,  and 90% credible interval of the estimates
lions_output_s.af_estimates = data.frame(parameter = colnames(lions_output_s.af), 
                                         mean = apply(lions_output_s.af, 2, mean),
                                         median = apply(lions_output_s.af, 2, median),
                                         lowerCI = apply(lions_output_s.af, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                         upperCI = apply(lions_output_s.af, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_s.af_estimates$parameter_clean = c("Base estimate (wet season in the grassland)",
                                                "Base estimate (dry season in the grassland)",
                                                "Age (dry season)",
                                                "Age (wet season)",
                                                "Habitat (grassland, wet season)",
                                                "Habitat (grassland, dry season)",
                                                "Habitat (woodland, wet season)",
                                                "Habitat (woodland, dry season)",
                                                "Number of females in the pride (wet season)",
                                                "Number of females in the pride (dry season)",
                                                "Number of nomadic coalitions in the home range (wet season)",
                                                "Number of nomadic coalitions in the home range (dry season)")


# Add season and parameter category
lions_output_s.af_estimates$season = rep(c("Wet", "Dry"), 6)
lions_output_s.af_estimates$parameter_plot = rep(c("Mean survival", 
                                                   "Age", 
                                                   "Habitat (grassland)", 
                                                   "Habitat (woodland)",
                                                   "Number of females\nin the pride", 
                                                   "Number of nomadic\ncoalitions in the home range"), 
                                                 each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_s.af_estimates$parameter_plot))),
                               variable = rep(c("Mean survival", "Age", 
                                                "Habitat (grassland)", 
                                                "Habitat (woodland)", 
                                                "Number of females\nin the pride", 
                                                "Number of nomadic\ncoalitions in the home range"), 
                                              each = 2 * n),
                               posterior = c(lions_output_s.af[1:nrow(lions_output_s.af), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean survival", 
                                               "Number of females\nin the pride", 
                                               "Number of nomadic\ncoalitions in the home range", 
                                               "Habitat (woodland)", 
                                               "Age")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.3.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_s.af_estimates$parameter_plot,
                     season = lions_output_s.af_estimates$season,
                     mean = lions_output_s.af_estimates$mean,
                     median = lions_output_s.af_estimates$median,
                     BCI90_upper = lions_output_s.af_estimates$upperCI,
                     BCI90_lower = lions_output_s.af_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean survival", 
                                     "Number of females\nin the pride", 
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Age", 
                                     "Habitat (woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "AF_Survival.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

afsurv_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, aes(x = median, y = variable, # Median
                                 color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("Yes", "No"),
                     values = c(1, 0.4)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 1), limits = c(-7.3, 3)) +
  theme_general() +
  ggtitle("Adult-female survival") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

afsurv_plot

dev.off()


## 3.3.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
af.survival.estimate = function(season = "wet",
                                habitat = "grassland",
                                age = 0,
                                nb.af.pride = 0,
                                nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.af", 
                                         colnames(lions_output_GLMM))][, 1]
    
    beta.age = lions_output_GLMM[, grep("s.af.beta.age", 
                                        colnames(lions_output_GLMM))][, 1]
    
    beta.nb.af.pride = lions_output_GLMM[, grep("s.af.beta.nb.af", 
                                                colnames(lions_output_GLMM))][, 1]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.af.beta.nb.nm", 
                                                  colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_s.af.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.af.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.af", 
                                         colnames(lions_output_GLMM))][, 2]
    
    beta.age = lions_output_GLMM[, grep("s.af.beta.age", 
                                        colnames(lions_output_GLMM))][, 2]
    
    beta.nb.af.pride = lions_output_GLMM[, grep("s.af.beta.nb.af", 
                                                colnames(lions_output_GLMM))][, 2]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.af.beta.nb.nm", 
                                                  colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_s.af.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = 0
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.af.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  pred =  mean.surv + beta.age * age + beta.nb.af.pride * nb.af.pride + 
          beta.nb.nm.coal.hr * nb.nm.coal.hr + beta.habitat + median.epsilons
  
  return(pred)
  
}


# Average survival per habitat and season

# Empty prediction dataframe
af.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                          habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
af.season.habitat.estimates$pred = plogis(apply(apply(af.season.habitat.estimates, 
                                                      1, 
                                                      FUN = function(x) af.survival.estimate(season = x[1], 
                                                                                             habitat = x[2])), 
                                                2, median))
af.season.habitat.estimates$lwr = plogis(apply(apply(af.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) af.survival.estimate(season = x[1], 
                                                                                            habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
af.season.habitat.estimates$upr = plogis(apply(apply(af.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) af.survival.estimate(season = x[1],
                                                                                            habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_AF_Survival_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(af.season.habitat.estimates, aes(x = season, y = pred, 
                                             colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Adult-female\nsurvival probability") +
  labs(title = "Adult-female survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average survival with age

# Empty prediction dataframe
af.age.estimates = expand.grid(season = c("wet", "dry"),
                               age = (seq(2, 15) - mean(age.unscaled, na.rm = T)) / (2 * sd(age.unscaled, na.rm = T)),
                               pred = NA)

# Fill in prediction and credible intervals
af.age.estimates$pred = plogis(apply(apply(af.age.estimates, 
                                           1, 
                                           FUN = function(x) af.survival.estimate(season = x[1],
                                                                                  age = as.numeric(x[2]))), 
                                     2, median))
af.age.estimates$lwr = plogis(apply(apply(af.age.estimates, 
                                          1, 
                                          FUN = function(x) af.survival.estimate(season = x[1],
                                                                                 age = as.numeric(x[2]))), 
                                    2, FUN = function(x) quantile(x, probs = 0.05)))
af.age.estimates$upr = plogis(apply(apply(af.age.estimates, 
                                          1, 
                                          FUN = function(x) af.survival.estimate(season = x[1],
                                                                                 age = as.numeric(x[2]))), 
                                    2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_AF_Survival_Age.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(af.age.estimates, aes(x = age, y = pred, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(af.age.estimates$age),
                     labels = seq(2, 15)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Age") +
  ylab("Adult-female survival probability") +
  labs(title = "Adult-female survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


# Average survival with number of adult females in the pride

# Empty prediction dataframe
af.nb.af.estimates = expand.grid(season = c("wet", "dry"),
                                 nb.af.pride = (seq(1, 12) - mean(nb.af.pride.unscaled, na.rm = T)) / (2 * sd(nb.af.pride.unscaled, na.rm = T)),
                                 pred = NA)

# Fill in prediction and credible intervals
af.nb.af.estimates$pred = plogis(apply(apply(af.nb.af.estimates, 
                                             1, 
                                             FUN = function(x) af.survival.estimate(season = x[1],
                                                                                    nb.af.pride = as.numeric(x[2]))),
                                       2, median))
af.nb.af.estimates$lwr = plogis(apply(apply(af.nb.af.estimates,
                                            1, 
                                            FUN = function(x) af.survival.estimate(season = x[1],
                                                                                   nb.af.pride = as.numeric(x[2]))),
                                      2, FUN = function(x) quantile(x, probs = 0.05)))
af.nb.af.estimates$upr = plogis(apply(apply(af.nb.af.estimates, 
                                            1, 
                                            FUN = function(x) af.survival.estimate(season = x[1],
                                                                                   nb.af.pride = as.numeric(x[2]))),
                                      2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_AF_Survival_NbAFPride.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(af.nb.af.estimates, aes(x = nb.af.pride, y = pred, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(af.nb.af.estimates$nb.af.pride),
                     labels = seq(1, 12)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Number of females\nin the pride") +
  ylab("Adult-female\nsurvival probability") +
  labs(title = "Adult-female survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


# Average survival with number of nomadic coalitions in the home range

# Empty prediction dataframe
af.nb.nm.coal.estimates = expand.grid(season = c("wet", "dry"),
                                      nb.nm.coal.hr = (seq(0, 6) - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)),
                                      pred = NA)

# Fill in prediction and credible intervals
af.nb.nm.coal.estimates$pred = plogis(apply(apply(af.nb.nm.coal.estimates,
                                                  1, 
                                                  FUN = function(x) af.survival.estimate(season = x[1],
                                                                                         nb.nm.coal.hr = as.numeric(x[2]))),
                                            2, median))
af.nb.nm.coal.estimates$lwr = plogis(apply(apply(af.nb.nm.coal.estimates, 
                                                 1, 
                                                 FUN = function(x) af.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]))), 
                                           2, FUN = function(x) quantile(x, probs = 0.05)))
af.nb.nm.coal.estimates$upr = plogis(apply(apply(af.nb.nm.coal.estimates, 
                                                 1, 
                                                 FUN = function(x) af.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]))),
                                           2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_AF_Survival_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(af.nb.nm.coal.estimates, aes(x = nb.nm.coal.hr, y = pred, 
                                    colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(af.nb.nm.coal.estimates$nb.nm.coal.hr),
                     labels = seq(0, 6)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Adult-female survival probability") +
  labs(title = "Adult-female survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


## 7.4. Young male survival ----
# -------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_s.ym = lions_output_GLMM[, grep("s.ym.", 
                                             colnames(lions_output_GLMM))]
lions_output_s.ym.epsilons = lions_output_s.ym[, grep("epsilon", 
                                                      colnames(lions_output_s.ym))]
lions_output_s.ym = lions_output_s.ym[, - grep("epsilon", 
                                               colnames(lions_output_s.ym))]
lions_output_s.ym = lions_output_s.ym[, - grep("sigma", 
                                               colnames(lions_output_s.ym))]


# Mean and 90% credible interval of the estimates
lions_output_s.ym_estimates = data.frame(parameter = colnames(lions_output_s.ym), 
                                         mean = apply(lions_output_s.ym, 2, mean),
                                         median = apply(lions_output_s.ym, 2, median),
                                         lowerCI = apply(lions_output_s.ym, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                         upperCI = apply(lions_output_s.ym, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_s.ym_estimates$parameter_clean = c("Base estimate (wet season in the grassland)",
                                                "Base estimate (dry season in the grassland)",
                                                "Habitat (grassland, wet season)",
                                                "Habitat (grassland, dry season)",
                                                "Habitat (woodland, wet season)",
                                                "Habitat (woodland, dry season)",
                                                "Number of females in the pride (wet season)",
                                                "Number of females in the pride (dry season)",
                                                "Number of nomadic coalitions in the home range (wet season)",
                                                "Number of nomadic coalitions in the home range (dry season)")


# Add season and parameter category
lions_output_s.ym_estimates$season = rep(c("Wet", "Dry"), 5)
lions_output_s.ym_estimates$parameter_plot = rep(c("Mean survival", 
                                                   "Habitat (grassland)", 
                                                   "Habitat\n(woodland)", 
                                                   "Number of females\nin the pride", 
                                                   "Number of nomadic\ncoalitions in the home range"), 
                                                 each = 2)

# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples 

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_s.ym_estimates$parameter_plot))),
                               variable = rep(c("Mean survival", 
                                                "Habitat (grassland)", 
                                                "Habitat\n(woodland)", 
                                                "Number of females\nin the pride", 
                                                "Number of nomadic\ncoalitions in the home range"), 
                                              each = 2 * n),
                               posterior = c(lions_output_s.ym[1:nrow(lions_output_s.ym), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean survival", 
                                               "Number of females\nin the pride", 
                                               "Number of nomadic\ncoalitions in the home range", 
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


# 7.4.1. Plotting the effect sizes ----
# --------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_s.ym_estimates$parameter_plot,
                     season = lions_output_s.ym_estimates$season,
                     mean = lions_output_s.ym_estimates$mean,
                     median = lions_output_s.ym_estimates$median,
                     BCI90_upper = lions_output_s.ym_estimates$upperCI,
                     BCI90_lower = lions_output_s.ym_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean survival", 
                                     "Number of females\nin the pride", 
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}

# Plot effect size
png(filename = "YM_Survival.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ymsurv_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping zero",
                     labels = c("Yes", "No"),
                     values = c(1, 0.4)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 2), limits = c(-3.9, 10)) +
  theme_general() +
  ggtitle("Young-male survival") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ymsurv_plot

dev.off()


## 7.4.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
ym.survival.estimate = function(season = "wet",
                                habitat = "grassland", 
                                nb.af.pride = 0,
                                nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.ym", colnames(lions_output_GLMM))][, 1]
    
    beta.nb.af.pride = lions_output_GLMM[, grep("s.ym.beta.nb.af.pride", 
                                                colnames(lions_output_GLMM))][, 1]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.ym.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_s.ym.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.ym.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.ym", colnames(lions_output_GLMM))][, 2]
    
    beta.nb.af.pride = lions_output_GLMM[, grep("s.ym.beta.nb.af.pride", 
                                                colnames(lions_output_GLMM))][, 2]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.ym.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_s.ym.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = 0
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.ym.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred =  mean.surv + beta.nb.af.pride * nb.af.pride + 
          beta.nb.nm.coal.hr * nb.nm.coal.hr + beta.habitat + median.epsilons
  
  
  return(pred)
  
}


# Average survival per habitat and season

# Empty prediction dataframe
ym.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                          habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
ym.season.habitat.estimates$pred = plogis(apply(apply(ym.season.habitat.estimates,
                                                      1, 
                                                      FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                             habitat = x[2])), 
                                                2, median))
ym.season.habitat.estimates$lwr = plogis(apply(apply(ym.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                            habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
ym.season.habitat.estimates$upr = plogis(apply(apply(ym.season.habitat.estimates,
                                                     1,
                                                     FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                            habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_YM_Survival_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ym.season.habitat.estimates, aes(x = season, y = pred, 
                                        colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Young-male survival\nprobability") +
  labs(title = "Young-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average survival with number of adult females in the pride

# Empty prediction dataframe
ym.nb.af.estimates = expand.grid(season = c("wet", "dry"),
                                 nb.af.pride = (seq(1, 12) - mean(nb.af.pride.unscaled, na.rm = T)) / (2 * sd(nb.af.pride.unscaled, na.rm = T)),
                                 pred = NA)

# Fill in prediction and credible intervals
ym.nb.af.estimates$pred = plogis(apply(apply(ym.nb.af.estimates, 
                                             1, 
                                             FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                    nb.af.pride = as.numeric(x[2]))), 
                                       2, median))
ym.nb.af.estimates$lwr = plogis(apply(apply(ym.nb.af.estimates, 
                                            1, 
                                            FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                   nb.af.pride = as.numeric(x[2]))), 
                                      2, FUN = function(x) quantile(x, probs = 0.05)))
ym.nb.af.estimates$upr = plogis(apply(apply(ym.nb.af.estimates, 
                                            1, 
                                            FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                   nb.af.pride = as.numeric(x[2]))), 
                                      2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_YM_Survival_NbAFPride.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ym.nb.af.estimates, aes(x = nb.af.pride, y = pred, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(ym.nb.af.estimates$nb.af.pride),
                     labels = seq(1, 12)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Number of females\nin the pride") +
  ylab("Young-male\nsurvival probability") +
  labs(title = "Young-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


# Average survival with number of nomadic coalitions in the home range

# Empty prediction dataframe
ym.nb.nm.coal.estimates = expand.grid(season = c("wet", "dry"),
                                      nb.nm.coal.hr = (seq(0, 6) - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)),
                                      pred = NA)

# Fill in prediction and credible intervals
ym.nb.nm.coal.estimates$pred = plogis(apply(apply(ym.nb.nm.coal.estimates, 
                                                  1, 
                                                  FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                         nb.nm.coal.hr = as.numeric(x[2]))),
                                            2, median))
ym.nb.nm.coal.estimates$lwr = plogis(apply(apply(ym.nb.nm.coal.estimates,
                                                 1, 
                                                 FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]))), 
                                           2, FUN = function(x) quantile(x, probs = 0.05)))
ym.nb.nm.coal.estimates$upr = plogis(apply(apply(ym.nb.nm.coal.estimates, 
                                                 1, 
                                                 FUN = function(x) ym.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]))),
                                           2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_YM_Survival_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ym.nb.nm.coal.estimates, aes(x = nb.nm.coal.hr, y = pred, 
                                    colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(ym.nb.nm.coal.estimates$nb.nm.coal.hr),
                     labels = seq(0, 6)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Young-male\nsurvival probability") +
  labs(title = "Young-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8),
        axis.title = element_text(colour = colPlot, size = 9),
        legend.title = element_text(colour = colPlot, size = 9),
        legend.text = element_text(colour = colPlot, size = 8),
        plot.title = element_text(colour = colPlot))

dev.off()  


## 7.5. Nomadic male survival ----
# ---------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_s.nm = lions_output_GLMM[, grep("s.nm.", 
                                             colnames(lions_output_GLMM))]
lions_output_s.nm.epsilons = lions_output_s.nm[, grep("epsilon", 
                                                      colnames(lions_output_s.nm))]
lions_output_s.nm = lions_output_s.nm[, - grep("epsilon", 
                                               colnames(lions_output_s.nm))]
lions_output_s.nm = lions_output_s.nm[, - grep("sigma", 
                                               colnames(lions_output_s.nm))]


# Mean, median, and 90% credible interval of the estimates
lions_output_s.nm_estimates = data.frame(parameter = colnames(lions_output_s.nm), 
                                         mean = apply(lions_output_s.nm, 2, mean),
                                         median = apply(lions_output_s.nm, 2, median),
                                         lowerCI = apply(lions_output_s.nm, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                         upperCI = apply(lions_output_s.nm, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_s.nm_estimates$parameter_clean = c("Base estimate (wet season in the grassland)",
                                                "Base estimate (dry season in the grassland)",
                                                "Coalition size (wet season)",
                                                "Coalition size (dry season)",
                                                "Habitat (grassland, wet season)",
                                                "Habitat (grassland, dry season)",
                                                "Habitat (woodland, wet season)",
                                                "Habitat (woodland, dry season)")

# Add season and parameter category
lions_output_s.nm_estimates$season = rep(c("Wet", "Dry"), 4)
lions_output_s.nm_estimates$parameter_plot = rep(c("Mean survival", 
                                                   "Coalition size", 
                                                   "Habitat (grassland)", 
                                                   "Habitat\n(woodland)"), 
                                                 each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_s.nm_estimates$parameter_plot))),
                               variable = rep(c("Mean survival", 
                                                "Coalition size", 
                                                "Habitat (grassland)", 
                                                "Habitat\n(woodland)"), 
                                              each = 2 * n),
                               posterior = c(lions_output_s.nm[1:nrow(lions_output_s.nm), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean survival", "Coalition size", 
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.5.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_s.nm_estimates$parameter_plot,
                     season = lions_output_s.nm_estimates$season,
                     mean = lions_output_s.nm_estimates$mean,
                     median = lions_output_s.nm_estimates$median,
                     BCI90_upper = lions_output_s.nm_estimates$upperCI,
                     BCI90_lower = lions_output_s.nm_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean survival", "Coalition size", 
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))

# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "NM_Survival.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

nmsurv_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), position = 
                   position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 2), limits = c(-5.3, 10)) +
  theme_general() +
  ggtitle("Nomadic-male survival") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

nmsurv_plot

dev.off()


## 7.5.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
nm.survival.estimate = function(season = "wet",
                                habitat = "grassland",
                                coal.size = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.nm", 
                                         colnames(lions_output_GLMM))][, 1]
    
    beta.coal.size = lions_output_GLMM[, grep("s.nm.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_s.nm.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.nm.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.nm", 
                                         colnames(lions_output_GLMM))][, 2]
    
    beta.coal.size = lions_output_GLMM[, grep("s.nm.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_s.nm.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.nm.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
    
    # Calculate vital-rate prediction
    pred =  mean.surv + beta.coal.size * coal.size + beta.habitat + median.epsilons
    
  }
  
  return(pred)
  
}


# Average survival per habitat and season

# Empty prediction dataframe
nm.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                               habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
nm.season.habitat.estimates$pred = plogis(apply(apply(nm.season.habitat.estimates, 
                                                      1, 
                                                      FUN = function(x) nm.survival.estimate(season = x[1],
                                                                                             habitat = x[2])), 
                                                2, median))
nm.season.habitat.estimates$lwr = plogis(apply(apply(nm.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) nm.survival.estimate(season = x[1],
                                                                                            habitat = x[2])),
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
nm.season.habitat.estimates$upr = plogis(apply(apply(nm.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) nm.survival.estimate(season = x[1],
                                                                                            habitat = x[2])),
                                               2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_NM_Survival_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(nm.season.habitat.estimates, aes(x = season, y = pred, 
                                        colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Nomadic-male\nsurvival probability") +
  labs(title = "Nomadic-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average survival with coalition size

# Empty prediction dataframe
nm.coal.size.estimates = expand.grid(season = c("wet", "dry"),
                                     coal.size = (seq(1, 4) - mean(coal.size.unscaled, na.rm = T)) / (2 * sd(coal.size.unscaled, na.rm = T)),
                                      pred = NA)

# Fill in prediction and credible intervals
nm.coal.size.estimates$pred = plogis(apply(apply(nm.coal.size.estimates, 
                                                 1, 
                                                 FUN = function(x) nm.survival.estimate(season = x[1],
                                                                                        coal.size = as.numeric(x[2]))), 
                                           2, mean))
nm.coal.size.estimates$lwr = plogis(apply(apply(nm.coal.size.estimates, 
                                                1, 
                                                FUN = function(x) nm.survival.estimate(season = x[1],
                                                                                       coal.size = as.numeric(x[2]))),
                                          2, FUN = function(x) quantile(x, probs = 0.05)))
nm.coal.size.estimates$upr = plogis(apply(apply(nm.coal.size.estimates, 
                                                1, 
                                                FUN = function(x) nm.survival.estimate(season = x[1],
                                                                                       coal.size = as.numeric(x[2]))), 
                                          2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_NM_Survival_CoalSize.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(nm.coal.size.estimates, aes(x = coal.size, y = pred, 
                                   colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(nm.coal.size.estimates$coal.size),
                     labels = seq(1, 4)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Coalition size") +
  ylab("Nomadic-male\nsurvival probability") +
  labs(title = "Nomadic-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


## 7.6. Resident male survival ----
# ----------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_s.rm = lions_output_GLMM[, grep("s.rm.", 
                                             colnames(lions_output_GLMM))]
lions_output_s.rm.epsilons = lions_output_s.rm[, grep("epsilon", 
                                                      colnames(lions_output_s.rm))]
lions_output_s.rm = lions_output_s.rm[, - grep("epsilon", 
                                               colnames(lions_output_s.rm))]
lions_output_s.rm = lions_output_s.rm[, - grep("sigma", 
                                               colnames(lions_output_s.rm))]


# Mean, median, and 90% credible interval of the estimates
lions_output_s.rm_estimates = data.frame(parameter = colnames(lions_output_s.rm), 
                                         mean = apply(lions_output_s.rm, 2, mean),
                                         median = apply(lions_output_s.rm, 2, median),
                                         lowerCI = apply(lions_output_s.rm, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                         upperCI = apply(lions_output_s.rm, 2, 
                                                         FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_s.rm_estimates$parameter_clean = c("Mean survival\n(wet season in the grassland)",
                                                "Mean survival\n(dry season in the grassland)",
                                                "Coalition size\n(wet season)",
                                                "Coalition size\n(dry season)",
                                                "Habitat\n(grassland, wet season)",
                                                "Habitat\n(grassland, dry season)",
                                                "Habitat\n(woodland, wet season)",
                                                "Habitat\n(woodland, dry season)",
                                                "Number of nomadic\ncoalitions in the home range (wet season)",
                                                "Number of nomadic\ncoalitions in the home range (dry season)")


# Add season and parameter category
lions_output_s.rm_estimates$season = rep(c("Wet", "Dry"), 5)
lions_output_s.rm_estimates$parameter_plot = rep(c("Mean survival", 
                                                   "Coalition size", 
                                                   "Habitat (grassland)", 
                                                   "Habitat (woodland)", 
                                                   "Number of nomadic\ncoalitions in the home range"), 
                                                 each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_s.rm_estimates$parameter_plot))),
                               variable = rep(c("Mean survival", 
                                                "Coalition size", 
                                                "Habitat (grassland)",
                                                "Habitat (woodland)", 
                                                "Number of nomadic\ncoalitions in the home range"), 
                                              each = 2 * n),
                               posterior = c(lions_output_s.rm[1:nrow(lions_output_s.rm), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean survival", "Coalition size", 
                                               "Number of nomadic\ncoalitions in the home range", 
                                               "Habitat (woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.6.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_s.rm_estimates$parameter_plot,
                     season = lions_output_s.rm_estimates$season,
                     mean = lions_output_s.rm_estimates$mean,
                     median = lions_output_s.rm_estimates$median,
                     BCI90_upper = lions_output_s.rm_estimates$upperCI,
                     BCI90_lower = lions_output_s.rm_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean survival", "Coalition size", 
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Habitat (woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Prediction plot
png(filename = "RM_Survival.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

rmsurv_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot 
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 linewidth = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, aes(x = median, y = variable, # Median
                                 color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 8, 1), limits = c(-2, 4.3)) +
  theme_general() +
  ggtitle("Resident-male survival") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

rmsurv_plot

dev.off()


## 7.6.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
rm.survival.estimate = function(season = "wet",
                                habitat = "grassland",
                                coal.size = 0,
                                nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.rm", 
                                    colnames(lions_output_GLMM))][, 1]
    
    beta.coal.size = lions_output_GLMM[, grep("s.rm.beta.coal.size", 
                                         colnames(lions_output_GLMM))][, 1]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.rm.beta.nb.nm.coal.hr", 
                                             colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_s.rm.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.rm.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.surv = lions_output_GLMM[, grep("mu.s.rm", 
                                         colnames(lions_output_GLMM))][, 2]
    
    beta.coal.size = lions_output_GLMM[, grep("s.rm.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 2]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("s.rm.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_s.rm.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("s.rm.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
    
    # Calculate vital-rate prediction
    pred =  mean.surv + beta.coal.size * coal.size + 
            beta.nb.nm.coal.hr * coal.size + beta.habitat + median.epsilons
    
  }
  
  return(pred)
  
}


# Average survival per habitat and season

# Empty prediction dataframe
rm.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                          habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
rm.season.habitat.estimates$pred = plogis(apply(apply(rm.season.habitat.estimates, 
                                                      1, 
                                                      FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                             habitat = x[2])), 
                                                2, median))
rm.season.habitat.estimates$lwr = plogis(apply(apply(rm.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                            habitat = x[2])), 
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
rm.season.habitat.estimates$upr = plogis(apply(apply(rm.season.habitat.estimates, 
                                                     1, 
                                                     FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                            habitat = x[2])),
                                               2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_RM_Survival_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rm.season.habitat.estimates, aes(x = season, y = pred, 
                                             colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Resident-male\nsurvival probability") +
  labs(title = "Resident-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average survival with coalition size

# Empty prediction dataframe
rm.coal.size.estimates = expand.grid(season = c("wet", "dry"),
                                          coal.size = (seq(1, 4) - mean(coal.size.unscaled, na.rm = T)) / (2 * sd(coal.size.unscaled, na.rm = T)),
                                          pred = NA)

# Fill in prediction and credible intervals
rm.coal.size.estimates$pred = plogis(apply(apply(rm.coal.size.estimates, 
                                                 1, 
                                                 FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                        coal.size = as.numeric(x[2]))), 
                                           2, median))
rm.coal.size.estimates$lwr = plogis(apply(apply(rm.coal.size.estimates, 
                                                1, 
                                                FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                       coal.size = as.numeric(x[2]))), 
                                          2, FUN = function(x) quantile(x, probs = 0.05)))
rm.coal.size.estimates$upr = plogis(apply(apply(rm.coal.size.estimates,
                                                1, 
                                                FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                       coal.size = as.numeric(x[2]))),
                                          2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_RM_Survival_CoalSize.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rm.coal.size.estimates, aes(x = coal.size, y = pred, 
                                   colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(rm.coal.size.estimates$coal.size),
                     labels = seq(1, 4)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Coalition size\n ") +
  ylab("Resident-male\nsurvival probability") +
  labs(title = "Resident-male survival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


# Average survival with number of nomadic coalitions in the home range

# Empty prediction dataframe
rm.nb.nm.coal.estimates = expand.grid(season = c("wet", "dry"),
                                      nb.nm.coal.hr = (seq(0, 6) - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)),
                                      pred = NA)

# Fill in prediction and credible intervals
rm.nb.nm.coal.estimates$pred = plogis(apply(apply(rm.nb.nm.coal.estimates, 
                                                  1, 
                                                  FUN = function(x) rm.survival.estimate(season = x[1], 
                                                                                         nb.nm.coal.hr = as.numeric(x[2]))), 
                                            2, median))
rm.nb.nm.coal.estimates$lwr = plogis(apply(apply(rm.nb.nm.coal.estimates, 
                                                 1, 
                                                 FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]))),
                                           2, FUN = function(x) quantile(x, probs = 0.05)))
rm.nb.nm.coal.estimates$upr = plogis(apply(apply(rm.nb.nm.coal.estimates,
                                                 1, 
                                                 FUN = function(x) rm.survival.estimate(season = x[1],
                                                                                        nb.nm.coal.hr = as.numeric(x[2]))),
                                           2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_RM_Survival_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rm.nb.nm.coal.estimates, aes(x = nb.nm.coal.hr, y = pred, 
                                    colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(rm.nb.nm.coal.estimates$nb.nm.coal.hr),
                     labels = seq(0, 6)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Resident-male\nsurvival probability") +
  labs(title = "Resident-male\nsurvival") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()  


## 7.7. Young male emigration ----
# ---------------------------

# Subset model output data and remove epsilons and sigma
lions_output_emig.ym = lions_output_GLMM[, grep("emig.ym", 
                                                colnames(lions_output_GLMM))]
lions_output_emig.ym.epsilons = lions_output_emig.ym[, grep("epsilon", 
                                                            colnames(lions_output_emig.ym))]
lions_output_emig.ym = lions_output_emig.ym[, - grep("epsilon", 
                                                     colnames(lions_output_emig.ym))]
lions_output_emig.ym = lions_output_emig.ym[, - grep("sigma", 
                                                     colnames(lions_output_emig.ym))]


# Mean, median, and 90% credible interval of the estimates
lions_output_emig.ym_estimates = data.frame(parameter = colnames(lions_output_emig.ym), 
                                            mean = apply(lions_output_emig.ym, 
                                                         2, mean),
                                            median = apply(lions_output_emig.ym, 
                                                         2, median),
                                            lowerCI = apply(lions_output_emig.ym, 
                                                            2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                            upperCI = apply(lions_output_emig.ym, 
                                                            2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_emig.ym_estimates$parameter_clean = c("Habitat (grassland, wet season)",
                                                   "Habitat (grassland, dry season)",
                                                   "Habitat (woodland, wet season)",
                                                   "Habitat (woodland, dry season)",
                                                   "Base estimate (wet season in the grassland)",
                                                   "Base estimate (dry season in the grassland)")

# Add season and parameter category
lions_output_emig.ym_estimates$season = rep(c("Wet", "Dry"), 3)
lions_output_emig.ym_estimates$parameter_plot = rep(c("Habitat (grassland)", "Habitat\n(woodland)", "Mean emigration\nprobability"), each = 2)


# Get posterior distribution for density plots
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_emig.ym_estimates$parameter_plot))),
                               variable = rep(c("Habitat (grassland)", 
                                                "Habitat\n(woodland)", 
                                                "Mean emigration\nprobability"), 
                                              each = 2 * n),
                               posterior = c(lions_output_emig.ym[1:nrow(lions_output_emig.ym), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean emigration\nprobability", 
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))



## 7.7.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_emig.ym_estimates$parameter_plot,
                     season = lions_output_emig.ym_estimates$season,
                     mean = lions_output_emig.ym_estimates$mean,
                     median = lions_output_emig.ym_estimates$median,
                     BCI90_upper = lions_output_emig.ym_estimates$upperCI,
                     BCI90_lower = lions_output_emig.ym_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean emigration\nprobability", 
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Prediction plot
png(filename = "YM_Emigration.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ymemig_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line 
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 1), limits = c(-4.3, 1.7)) +
  theme_general() +
  ggtitle("Young-male emigration") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ymemig_plot

dev.off()


## 7.7.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
ym.emig.estimate = function(season = "wet",
                            habitat = "grassland"){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.emig = lions_output_GLMM[, grep("mu.emig.ym", 
                                         colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_emig.ym.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("emig.ym.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.emig = lions_output_GLMM[, grep("mu.emig.ym", 
                                         colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_emig.ym.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("emig.ym.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred =  mean.emig + beta.habitat + median.epsilons
  
  return(pred)
  
}


# Average emigration probability per habitat and season

# Empty prediction dataframe
ym.emig.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                               habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
ym.emig.season.habitat.estimates$pred = plogis(apply(apply(ym.emig.season.habitat.estimates, 
                                                           1, 
                                                           FUN = function(x) ym.emig.estimate(season = x[1],
                                                                                              habitat = x[2])), 
                                                     2, median))
ym.emig.season.habitat.estimates$lwr = plogis(apply(apply(ym.emig.season.habitat.estimates, 
                                                          1,
                                                          FUN = function(x) ym.emig.estimate(season = x[1],
                                                                                             habitat = x[2])), 
                                                    2, FUN = function(x) quantile(x, probs = 0.05)))
ym.emig.season.habitat.estimates$upr = plogis(apply(apply(ym.emig.season.habitat.estimates, 
                                                          1, 
                                                          FUN = function(x) ym.emig.estimate(season = x[1],
                                                                                             habitat = x[2])), 
                                                    2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_YM_Emigration_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ym.emig.season.habitat.estimates, aes(x = season, y = pred, 
                                             colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Young-male\nemigration probability") +
  labs(title = "Young-male emigration") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


## 7.8. Young male transition to nomadic male (once emigrated) ----
# ------------------------------------------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_t.ym.nm = lions_output_GLMM[, grep("t.ym.nm.", 
                                                colnames(lions_output_GLMM))]
lions_output_t.ym.nm.epsilons = lions_output_t.ym.nm[, grep("epsilon", 
                                                            colnames(lions_output_t.ym.nm))]
lions_output_t.ym.nm = lions_output_t.ym.nm[, - grep("epsilon", 
                                                     colnames(lions_output_t.ym.nm))]
lions_output_t.ym.nm = lions_output_t.ym.nm[, - grep("sigma", 
                                                     colnames(lions_output_t.ym.nm))]


# Mean, median, and 90% credible interval of the estimates
lions_output_t.ym.nm_estimates = data.frame(parameter = colnames(lions_output_t.ym.nm), 
                                            mean = apply(lions_output_t.ym.nm, 
                                                         2, mean),
                                            median = apply(lions_output_t.ym.nm, 
                                                         2, median),
                                            lowerCI = apply(lions_output_t.ym.nm, 
                                                            2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                            upperCI = apply(lions_output_t.ym.nm, 
                                                            2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_t.ym.nm_estimates$parameter_clean = c("Base estimate (wet season in the grassland)",
                                                   "Base estimate (dry season in the grassland)")


# Add season and parameter category
lions_output_t.ym.nm_estimates$season = rep(c("Wet", "Dry"), 1)
lions_output_t.ym.nm_estimates$parameter_plot = rep(c("Mean transition\nprobability"), each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_t.ym.nm_estimates$parameter_plot))),
                               variable = rep(c("Mean transition\nprobability"), 
                                              each = 2 * n),
                               posterior = c(lions_output_t.ym.nm[1:nrow(lions_output_t.ym.nm), ]))

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean transition\nprobability"))
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.8.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_t.ym.nm_estimates$parameter_plot,
                     season = lions_output_t.ym.nm_estimates$season,
                     mean = lions_output_t.ym.nm_estimates$mean,
                     median = lions_output_t.ym.nm_estimates$median,
                     BCI90_upper = lions_output_t.ym.nm_estimates$upperCI,
                     BCI90_lower = lions_output_t.ym.nm_estimates$lowerCI,
                     overlapping0 = FALSE)

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable,
                          levels = unique(df.plot$variable)[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "YM_TransitionToNM.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ymtrans_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season),
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 2), limits = c(-2.5, 10)) +
  theme_general() +
  ggtitle("Young-male transition to nomadic male (once emigrated)") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ymtrans_plot

dev.off()


## 7.8.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
ym.trans.estimate = function(season = "wet"){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.trans = lions_output_GLMM[, grep("mu.t.ym.nm", 
                                        colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_t.ym.nm.epsilons[, seq(1, 59, 2)], 1, median)
    
  }
  
  else{
    
    # Get dry-season parameter values
    mean.trans = lions_output_GLMM[, grep("mu.t.ym.nm", 
                                          colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons.dry = apply(lions_output_t.ym.nm.epsilons[, seq(2, 60, 2)], 1, median)
    
  }
  
  # Calculate vital-rate prediction
  pred = mean.trans + median.epsilons
  
  return(pred)
  
}


# Average transition probability per habitat and season

# Empty prediction dataframe
ym.trans.season.habitat.estimates = expand.grid(season = c("wet", "dry"))


# Fill in prediction and credible intervals
ym.trans.season.habitat.estimates$pred = plogis(apply(apply(ym.trans.season.habitat.estimates, 
                                                            1, 
                                                            FUN = function(x) ym.trans.estimate(season = x[1])),
                                                      2, median))
ym.trans.season.habitat.estimates$lwr = plogis(apply(apply(ym.trans.season.habitat.estimates, 
                                                           1, 
                                                           FUN = function(x) ym.trans.estimate(season = x[1])),
                                                     2, FUN = function(x) quantile(x, probs = 0.05)))
ym.trans.season.habitat.estimates$upr = plogis(apply(apply(ym.trans.season.habitat.estimates, 
                                                           1, 
                                                           FUN = function(x) ym.trans.estimate(season = x[1])), 
                                                     2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_YM_TransitionNM_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(ym.trans.mean.season.habitat.estimates, aes(x = season, y = pred)) +
  geom_point(position = position_dodge(0.8), 
             colour = colPlot) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5, colour = colPlot) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  xlab("Season") +
  ylab("Young-male to nomad\ntransition probability") +
  labs(title = "Young-male transition to nomad") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


## 7.9. Nomadic male takeover ----
# ----------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_takeover = lions_output_GLMM[, grep("takeover.", 
                                                 colnames(lions_output_GLMM))]
lions_output_takeover.epsilons = lions_output_takeover[, grep("epsilon",
                                                              colnames(lions_output_takeover))]
lions_output_takeover = lions_output_takeover[, - grep("epsilon", 
                                                       colnames(lions_output_takeover))]
lions_output_takeover = lions_output_takeover[, - grep("sigma", 
                                                       colnames(lions_output_takeover))]


# Mean, median, and 90% credible interval of the estimates
lions_output_takeover_estimates = data.frame(parameter = colnames(lions_output_takeover), 
                                             mean = apply(lions_output_takeover, 2, mean),
                                             median = apply(lions_output_takeover, 2, median),
                                             lowerCI = apply(lions_output_takeover, 2, 
                                                             FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                             upperCI = apply(lions_output_takeover, 2, 
                                                             FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_takeover_estimates$parameter_clean = c("Base estimate (wet season in the grassland)",
                                                    "Base estimate (dry season in the grassland)",
                                                    "Coalition size (wet season)",
                                                    "Coalition size (dry season)",
                                                    "Habitat (grassland, wet season)",
                                                    "Habitat (grassland, dry season)",
                                                    "Habitat (woodland, wet season)",
                                                    "Habitat (woodland, dry season)")


# Add season and parameter category
lions_output_takeover_estimates$season = rep(c("Wet", "Dry"), 4)
lions_output_takeover_estimates$parameter_plot = rep(c("Mean takeover\nprobability", 
                                                       "Coalition size", 
                                                       "Habitat (grassland)", 
                                                       "Habitat\n(woodland)"), 
                                                     each = 2)


# Get posterior distribution for density
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_takeover_estimates$parameter_plot))),
                               variable = rep(c("Mean takeover\nprobability", 
                                                "Coalition size", 
                                                "Habitat (grassland)",
                                                "Habitat\n(woodland)"),
                                              each = 2 * n),
                               posterior = c(lions_output_takeover[1:nrow(lions_output_takeover), ]))

# Remove 0s to discard the habitat effect for the plains, which is the baseline
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean takeover\nprobability",
                                               "Coalition size", 
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season,
                                  levels = unique(df.plot.posterior$season))


## 7.9.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_takeover_estimates$parameter_plot,
                     season = lions_output_takeover_estimates$season,
                     mean = lions_output_takeover_estimates$mean,
                     median = lions_output_takeover_estimates$median,
                     BCI90_upper = lions_output_takeover_estimates$upperCI,
                     BCI90_lower = lions_output_takeover_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean takeover\nprobability", 
                                     "Coalition size",
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "NM_Takeover.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

nmtakeover_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, aes(x = median, y = variable, # Median
                                 color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 1), limits = c(-3.7, 4)) +
  theme_general() +
  ggtitle("Nomadic-male takeover") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

nmtakeover_plot

dev.off()


## 7.9.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
nm.takeover.estimate = function(season = "wet",
                                habitat = "grassland",
                                coal.size = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.takeover = lions_output_GLMM[, grep("mu.takeover", 
                                             colnames(lions_output_GLMM))][, 1]
    
    beta.coal.size = lions_output_GLMM[, grep("takeover.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_takeover.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat.wet = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("takeover.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.takeover = lions_output_GLMM[, grep("mu.takeover", 
                                             colnames(lions_output_GLMM))][, 2]
    
    beta.coal.size = lions_output_GLMM[, grep("takeover.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 2]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_takeover.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("takeover.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate prediction
  pred =  mean.takeover + beta.coal.size * coal.size + 
          beta.habitat + median.epsilons
  
  return(pred)
  
}


# Average takeover probability per habitat and season

# Empty prediction dataframe
nm.takeover.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                                   habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals 
nm.takeover.season.habitat.estimates$pred = plogis(apply(apply(nm.takeover.season.habitat.estimates,
                                                               1, 
                                                               FUN = function(x) nm.takeover.estimate(season = x[1],
                                                                                                      habitat = x[2])), 
                                                         2, median))
nm.takeover.season.habitat.estimates$lwr = plogis(apply(apply(nm.takeover.season.habitat.estimates, 
                                                              1, 
                                                              FUN = function(x) nm.takeover.estimate(season = x[1], 
                                                                                                     habitat = x[2])),
                                                        2, FUN = function(x) quantile(x, probs = 0.05)))
nm.takeover.season.habitat.estimates$upr = plogis(apply(apply(nm.takeover.season.habitat.estimates, 
                                                              1, 
                                                              FUN = function(x) nm.takeover.estimate(season = x[1],
                                                                                                     habitat = x[2])), 
                                                        2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_NM_Takeover_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(nm.takeover.season.habitat.estimates, aes(x = season, y = pred,
                                                 colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Nomadic-male\ntakeover probability") +
  labs(title = "Nomadic-male takeover") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average takeover probability with coalition size

# Empty prediction dataframe
nm.takeover.coal.size.estimates = expand.grid(season = c("wet", "dry"),
                                              coal.size = (seq(1, 4) - mean(coal.size.unscaled, na.rm = T)) / (2 * sd(coal.size.unscaled, na.rm = T)),
                                              pred = NA)

# Fill in prediction and credible intervals
nm.takeover.coal.size.estimates$pred = plogis(apply(apply(nm.takeover.coal.size.estimates, 
                                                          1, 
                                                          FUN = function(x) nm.takeover.estimate(season = x[1],
                                                                                                 coal.size = as.numeric(x[2]))), 
                                                    2, median))
nm.takeover.coal.size.estimates$lwr = plogis(apply(apply(nm.takeover.coal.size.estimates, 
                                                         1, 
                                                         FUN = function(x) nm.takeover.estimate(season = x[1],
                                                                                                coal.size = as.numeric(x[2]))), 
                                                   2, FUN = function(x) quantile(x, probs = 0.05)))
nm.takeover.coal.size.estimates$upr = plogis(apply(apply(nm.takeover.coal.size.estimates, 
                                                         1,
                                                         FUN = function(x) nm.takeover.estimate(season = x[1],
                                                                                                coal.size = as.numeric(x[2]))), 
                                                   2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_NM_Takeover_CoalSize.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(nm.takeover.coal.size.estimates, aes(x = coal.size, y = pred, 
                                            colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(nm.takeover.coal.size.estimates$coal.size),
                     labels = seq(1, 4)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Coalition size\n ") +
  ylab("Nomadic-male\ntakeover probability") +
  labs(title = "Nomadic-male takeover") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


## 7.10. Resident male eviction ----
# -----------------------------

# Subset model output data and remove epsilons
lions_output_eviction = lions_output_GLMM[, grep("eviction",
                                                 colnames(lions_output_GLMM))]
lions_output_eviction.epsilons = lions_output_eviction[, grep("epsilon", 
                                                              colnames(lions_output_eviction))]
lions_output_eviction = lions_output_eviction[, - grep("epsilon", 
                                                       colnames(lions_output_eviction))]
lions_output_eviction = lions_output_eviction[, - grep("sigma",
                                                       colnames(lions_output_eviction))]


# Mean, median, and 90% credible interval of the estimates
lions_output_eviction_estimates = data.frame(parameter = colnames(lions_output_eviction), 
                                             mean = apply(lions_output_eviction, 
                                                          2, mean),
                                             median = apply(lions_output_eviction, 
                                                          2, median),
                                             lowerCI = apply(lions_output_eviction, 
                                                             2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                             upperCI = apply(lions_output_eviction, 
                                                             2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_eviction_estimates$parameter_clean = c("Coalition size (wet season)",
                                                    "Coalition size (dry season)",
                                                    "Habitat (grassland, wet season)",
                                                    "Habitat (grassland, dry season)",
                                                    "Habitat (woodland, wet season)",
                                                    "Habitat (woodland, dry season)",
                                                    "Number of nomadic coalitions in the home range (wet season)",
                                                    "Number of nomadic coalitions in the home range (dry season)",
                                                    "Base estimate (wet season in the grassland)",
                                                    "Base estimate (dry season in the grassland)")


# Add season and parameter category
lions_output_eviction_estimates$season = rep(c("Wet", "Dry"), 5)
lions_output_eviction_estimates$parameter_plot = rep(c("Coalition size", 
                                                       "Habitat (grassland)", 
                                                       "Habitat\n(woodland)", 
                                                       "Number of nomadic\ncoalitions in the home range",
                                                       "Mean eviction\nprobability"), 
                                                     each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_eviction_estimates$parameter_plot))),
                               variable = rep(c("Coalition size", 
                                                "Habitat (grassland)",
                                                "Habitat\n(woodland)", 
                                                "Number of nomadic\ncoalitions in the home range",
                                                "Mean eviction\nprobability"), 
                                              each = 2 * n),
                               posterior = c(lions_output_eviction[1:nrow(lions_output_eviction), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean eviction\nprobability",
                                               "Coalition size", 
                                               "Number of nomadic\ncoalitions in the home range",
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season,
                                  levels = unique(df.plot.posterior$season))


## 7.9.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_eviction_estimates$parameter_plot,
                     season = lions_output_eviction_estimates$season,
                     mean = lions_output_eviction_estimates$mean,
                     median = lions_output_eviction_estimates$median,
                     BCI90_upper = lions_output_eviction_estimates$upperCI,
                     BCI90_lower = lions_output_eviction_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean eviction\nprobability", 
                                     "Coalition size", 
                                     "Number of nomadic\ncoalitions in the home range",
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "RM_Eviction.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

rmeviction_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plots
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3,
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals 
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7),
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, aes(x = median, y = variable, # Median
                                 color = season, alpha = overlapping0),
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 2), limits = c(-10, 5.6)) +
  theme_general() +
  ggtitle("Resident-male eviction") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

rmeviction_plot

dev.off()


## 7.9.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
rm.eviction.estimate = function(season = "wet",
                                habitat = "grassland",
                                coal.size = 0,
                                nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.eviction = lions_output_GLMM[, grep("mu.eviction", 
                                             colnames(lions_output_GLMM))][, 1]
    
    beta.coal.size = lions_output_GLMM[, grep("eviction.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 1]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("eviction.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_eviction.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("eviction.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.eviction = lions_output_GLMM[, grep("mu.eviction", 
                                             colnames(lions_output_GLMM))][, 2]
    
    beta.coal.size = lions_output_GLMM[, grep("eviction.beta.coal.size", 
                                              colnames(lions_output_GLMM))][, 2]
    
    beta.nb.nm.coal.hr = lions_output_GLMM[, grep("eviction.beta.nb.nm.coal.hr", 
                                                  colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_eviction.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("eviction.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred =  mean.eviction + beta.coal.size * coal.size + 
          beta.nb.nm.coal.hr * nb.nm.coal.hr + beta.habitat + median.epsilons
  
  return(pred)
  
}


# Average eviction probability per habitat and season

# Empty prediction dataframe
rm.eviction.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                                        habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
rm.eviction.season.habitat.estimates$pred = plogis(apply(apply(rm.eviction.season.habitat.estimates, 
                                                               1,
                                                               FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                      habitat = x[2])),
                                                         2, median))
rm.eviction.season.habitat.estimates$lwr = plogis(apply(apply(rm.eviction.season.habitat.estimates, 
                                                              1, 
                                                              FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                     habitat = x[2])), 
                                                        2, FUN = function(x) quantile(x, probs = 0.05)))
rm.eviction.season.habitat.estimates$upr = plogis(apply(apply(rm.eviction.season.habitat.estimates, 
                                                              1, 
                                                              FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                     habitat = x[2])),
                                                        2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_RM_Eviction_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rm.eviction.season.habitat.estimates, aes(x = season, y = pred, 
                                                 colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Resident-male\neviction probability") +
  labs(title = "Resident-male eviction") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average eviction probability with coalition size

# Empty prediction dataframe
rm.eviction.coal.size.estimates = expand.grid(season = c("wet", "dry"),
                                                   coal.size = (seq(1, 4) - mean(coal.size.unscaled, na.rm = T)) / (2 * sd(coal.size.unscaled, na.rm = T)),
                                                   pred = NA)

# Fill in prediction and credible intervals
rm.eviction.coal.size.estimates$pred = plogis(apply(apply(rm.eviction.coal.size.estimates,
                                                          1,
                                                          FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                 coal.size = as.numeric(x[2]))), 
                                                    2, median))
rm.eviction.coal.size.estimates$lwr = plogis(apply(apply(rm.eviction.coal.size.estimates, 
                                                         1, 
                                                         FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                coal.size = as.numeric(x[2]))),
                                                   2, FUN = function(x) quantile(x, probs = 0.05)))
rm.eviction.coal.size.estimates$upr = plogis(apply(apply(rm.eviction.coal.size.estimates, 
                                                         1, 
                                                         FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                coal.size = as.numeric(x[2]))), 
                                                   2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_RM_Eviction_CoalSize.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rm.eviction.coal.size.estimates, aes(x = coal.size, y = pred, 
                                            colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(rm.eviction.coal.size.estimates$coal.size),
                     labels = seq(1, 4)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Coalition size") +
  ylab("Resident-male\neviction probability") +
  labs(title = "Resident-male eviction") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


# Average eviction probability with number of nomadic coalitions in the home range

# Empty prediction dataframe
rm.eviction.nb.nm.coal.estimates = expand.grid(season = c("wet", "dry"),
                                                    nm.nm.coal.hr = (seq(0, 6) - mean(coal.size.unscaled, na.rm = T)) / (2 * sd(coal.size.unscaled, na.rm = T)),
                                                    pred = NA)

# Fill in prediction and credible intervals
rm.eviction.nb.nm.coal.estimates$pred = plogis(apply(apply(rm.eviction.nb.nm.coal.estimates, 
                                                           1, 
                                                           FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                  nb.nm.coal.hr = as.numeric(x[2]))),
                                                     2, median))
rm.eviction.nb.nm.coal.estimates$lwr = plogis(apply(apply(rm.eviction.nb.nm.coal.estimates,
                                                          1, 
                                                          FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                 nb.nm.coal.hr = as.numeric(x[2]))), 
                                                    2, FUN = function(x) quantile(x, probs = 0.05)))
rm.eviction.nb.nm.coal.estimates$upr = plogis(apply(apply(rm.eviction.nb.nm.coal.estimates, 
                                                          1,
                                                          FUN = function(x) rm.eviction.estimate(season = x[1],
                                                                                                 nb.nm.coal.hr = as.numeric(x[2]))),
                                                    2, FUN = function(x) quantile(x, probs = 0.95)))


png(filename = "Predictions_RM_Eviction_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rm.eviction.nb.nm.coal.estimates, aes(x = nm.nm.coal.hr, y = pred,
                                             colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, colour = season, fill = season), 
              alpha = 0.2) +
  scale_x_continuous(breaks = unique(rm.eviction.nb.nm.coal.estimates$nm.nm.coal.hr),
                     labels = seq(0, 6)) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Resident-male\neviction probability") +
  labs(title = "Resident-male eviction") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


## 7.10. Pride detection probability ----
# ----------------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_dp.pride = lions_output_GLMM[, grep("dp.pride", 
                                                 colnames(lions_output_GLMM))]
lions_output_dp.pride.epsilons = lions_output_dp.pride[, grep("epsilon", 
                                                              colnames(lions_output_dp.pride))]
lions_output_dp.pride = lions_output_dp.pride[, - grep("epsilon",
                                                       colnames(lions_output_dp.pride))]
lions_output_dp.pride = lions_output_dp.pride[, - grep("sigma",
                                                       colnames(lions_output_dp.pride))]


# Mean and 90% credible interval of estimates
lions_output_dp.pride_estimates = data.frame(parameter = colnames(lions_output_dp.pride), 
                                             mean = apply(lions_output_dp.pride, 
                                                          2, mean),
                                             median = apply(lions_output_dp.pride, 
                                                          2, median),
                                             lowerCI = apply(lions_output_dp.pride, 
                                                             2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                             upperCI = apply(lions_output_dp.pride, 
                                                             2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_dp.pride_estimates$parameter_clean = c("Habitat (grassland, wet season)",
                                                   "Habitat (grassland, dry season)",
                                                   "Habitat (woodland, wet season)",
                                                   "Habitat (woodland, dry season)",
                                                   "Base estimate (wet season in the grassland)",
                                                   "Base estimate (dry season in the grassland)")


# Add season and parameter category
lions_output_dp.pride_estimates$season = rep(c("Wet", "Dry"), 3)
lions_output_dp.pride_estimates$parameter_plot = rep(c("Habitat (grassland)",
                                                       "Habitat\n(woodland)", 
                                                       "Mean detection\nprobability"),
                                                     each = 2)


# Get posterior distribution for density plots
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_dp.pride_estimates$parameter_plot))),
                               variable = rep(c("Habitat (grassland)", 
                                                "Habitat\n(woodland)",
                                                "Mean detection\nprobability"),
                                              each = 2 * n),
                               posterior = c(lions_output_dp.pride[1:nrow(lions_output_dp.pride), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean detection\nprobability", 
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.10.1. Plotting the effect sizes ----
# ----------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_dp.pride_estimates$parameter_plot,
                     season = lions_output_dp.pride_estimates$season,
                     mean = lions_output_dp.pride_estimates$mean,
                     median = lions_output_dp.pride_estimates$median,
                     BCI90_upper = lions_output_dp.pride_estimates$upperCI,
                     BCI90_lower = lions_output_dp.pride_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean detection\nprobability", 
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season,
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "DP_Pride.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

dppride_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season, alpha = overlapping0), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 1), limits = c(-3.5, 5.4)) +
  theme_general() +
  ggtitle("Pride detection probability") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

dppride_plot

dev.off()


## 7.10.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
dp.pride.estimate = function(season = "wet",
                             habitat = "grassland"){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.dp = lions_output_GLMM[, grep("mu.dp.pride", 
                                       colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_dp.pride.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("dp.pride.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.dp = lions_output_GLMM[, grep("mu.dp.pride", 
                                       colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_dp.pride.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("dp.pride.beta.habitat.woodland", 
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred = mean.dp + beta.habitat + median.epsilons
  
  return(pred)
  
}


# Average detection probability per habitat and season

# Empty prediction dataframe
dp.pride.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                                     habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
dp.pride.season.habitat.estimates$pred = plogis(apply(apply(dp.pride.season.habitat.estimates, 
                                                            1, 
                                                            FUN = function(x) dp.pride.estimate(season = x[1], 
                                                                                                habitat = x[2])), 
                                                      2, median))
dp.pride.season.habitat.estimates$lwr = plogis(apply(apply(dp.pride.season.habitat.estimates,
                                                           1, 
                                                           FUN = function(x) dp.pride.estimate(season = x[1],
                                                                                               habitat = x[2])),
                                                     2, FUN = function(x) quantile(x, probs = 0.05)))
dp.pride.season.habitat.estimates$upr = plogis(apply(apply(dp.pride.season.habitat.estimates,
                                                           1,
                                                           FUN = function(x) dp.pride.estimate(season = x[1],
                                                                                               habitat = x[2])),
                                                     2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_DP_Pride_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(dp.pride.season.habitat.estimates, aes(x = season, y = pred, 
                                              colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Pride detection probability") +
  labs(title = "Pride detection probability") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


## 7.11. Nomad detection probability ----
# ----------------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_dp.nm = lions_output_GLMM[, grep("dp.nm", 
                                              colnames(lions_output_GLMM))]
lions_output_dp.nm.epsilons = lions_output_dp.nm[, grep("epsilon", 
                                                        colnames(lions_output_dp.nm))]

lions_output_dp.nm = lions_output_dp.nm[, - grep("epsilon", 
                                                 colnames(lions_output_dp.nm))]
lions_output_dp.nm = lions_output_dp.nm[, - grep("sigma", 
                                                 colnames(lions_output_dp.nm))]


# Mean, median, and 90% credible interval of the estimates
lions_output_dp.nm_estimates = data.frame(parameter = colnames(lions_output_dp.nm), 
                                          mean = apply(lions_output_dp.nm, 
                                                       2, mean),
                                          median = apply(lions_output_dp.nm, 
                                                       2, median),
                                          lowerCI = apply(lions_output_dp.nm, 
                                                          2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                          upperCI = apply(lions_output_dp.nm, 
                                                          2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_dp.nm_estimates$parameter_clean = c("Habitat (grassland, wet season)",
                                                 "Habitat (grassland, dry season)",
                                                 "Habitat (woodland, wet season)",
                                                 "Habitat (woodland, dry season)",
                                                 "Base estimate (wet season in the grassland)",
                                                 "Base estimate (dry season in the grassland)")


# Add season and parameter category
lions_output_dp.nm_estimates$season = rep(c("Wet", "Dry"), 3)
lions_output_dp.nm_estimates$parameter_plot = rep(c("Habitat (grassland)", 
                                                    "Habitat\n(woodland)",
                                                    "Mean detection\nprobability"),
                                                  each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_dp.nm_estimates$parameter_plot))),
                               variable = rep(c("Habitat (grassland)", 
                                                "Habitat\n(woodland)",
                                                "Mean detection\nprobability"), 
                                              each = 2 * n),
                               posterior = c(lions_output_dp.nm[1:nrow(lions_output_dp.nm), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean detection\nprobability",
                                               "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.11.1. Plotting the effect sizes ----
# -----------------å----------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_dp.nm_estimates$parameter_plot,
                     season = lions_output_dp.nm_estimates$season,
                     mean = lions_output_dp.nm_estimates$mean,
                     median = lions_output_dp.nm_estimates$median,
                     BCI90_upper = lions_output_dp.nm_estimates$upperCI,
                     BCI90_lower = lions_output_dp.nm_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean detection\nprobability", 
                                     "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Factorize parameter names and season
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effe size
png(filename = "DP_NM.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

dpnm_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season),
                       color = NA, size = 1, alpha = 0.3,
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season, alpha = overlapping0), 
                 position = position_dodge(width = 0.7), 
                 size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season, alpha = overlapping0),
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  scale_alpha_manual(name = "Overlapping 0",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 2), limits = c(-9.3, 10)) +
  theme_general() +
  ggtitle("Nomadic male detection probability") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

dpnm_plot

dev.off()


## 7.11.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
dp.nm.estimate = function(season = "wet",
                          habitat = "grassland"){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.dp = lions_output_GLMM[, grep("mu.dp.nm", 
                                       colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_dp.nm.epsilons[, seq(1, 59, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("dp.nm.beta.habitat.woodland",
                                              colnames(lions_output_GLMM))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.dp = lions_output_GLMM[, grep("mu.dp.nm", 
                                       colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_dp.nm.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_GLMM))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_GLMM[, grep("dp.nm.beta.habitat.woodland",
                                              colnames(lions_output_GLMM))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred =  mean.dp + beta.habitat + median.epsilons
  
  return(pred)
  
}


# Average detection probability per habitat and season

# Empty prediction dataframe
dp.nm.season.habitat.estimates = expand.grid(season = c("wet", "dry"),
                                             habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
dp.nm.season.habitat.estimates$pred = plogis(apply(apply(dp.nm.season.habitat.estimates, 
                                                         1, 
                                                         FUN = function(x) dp.nm.estimate(season = x[1],
                                                                                          habitat = x[2])),
                                                   2, median))
dp.nm.season.habitat.estimates$lwr = plogis(apply(apply(dp.nm.season.habitat.estimates,
                                                        1, 
                                                        FUN = function(x) dp.nm.estimate(season = x[1],
                                                                                         habitat = x[2])), 
                                                  2, FUN = function(x) quantile(x, probs = 0.05)))
dp.nm.season.habitat.estimates$upr = plogis(apply(apply(dp.nm.season.habitat.estimates,
                                                        1, 
                                                        FUN = function(x) dp.nm.estimate(season = x[1],
                                                                                         habitat = x[2])),
                                                  2, FUN = function(x) quantile(x, probs = 0.95)))

# Prediction plot
png(filename = "Predictions_DP_NM_Habitat_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(dp.nm.season.habitat.estimates, aes(x = season, y = pred, 
                                           colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Nomadic male\ndetection probability") +
  labs(title = "Nomadic male detection probability") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()


## 7.12. Death recovery probability ----
# ---------------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_dp.dead = lions_output_GLMM[, grep("dp.dead",
                                                colnames(lions_output_GLMM))]
lions_output_dp.dead.epsilons = lions_output_dp.dead[, grep("epsilon",
                                                            colnames(lions_output_dp.dead))]
lions_output_dp.dead = lions_output_dp.dead[, - grep("epsilon", 
                                                     colnames(lions_output_dp.dead))]
lions_output_dp.dead = lions_output_dp.dead[, - grep("sigma", 
                                                     colnames(lions_output_dp.dead))]


# Mean, median, and 90% credible interval of the estimates
lions_output_dp.dead_estimates = data.frame(parameter = colnames(lions_output_dp.dead), 
                                            mean = apply(lions_output_dp.dead, 
                                                         2, mean),
                                            median = apply(lions_output_dp.dead, 
                                                         2, median),
                                            lowerCI = apply(lions_output_dp.dead,
                                                            2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                            upperCI = apply(lions_output_dp.dead, 
                                                            2, FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_dp.dead_estimates$parameter_clean = c("Base estimate (wet season in the grassland)",
                                                   "Base estimate (dry season in the grassland)")


# Add season and parameter category
lions_output_dp.dead_estimates$season = rep(c("Wet", "Dry"), 1)
lions_output_dp.dead_estimates$parameter_plot = rep(c("Mean detection\nprobability"), 
                                                    each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_GLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n),
                                            length(unique(lions_output_dp.dead_estimates$parameter_plot))),
                               variable = rep(c("Mean detection\nprobability"),
                                              each = 2 * n),
                               posterior = c(lions_output_dp.dead[1:nrow(lions_output_dp.dead), ]))

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable,
                                    levels = c("Mean detection\nprobability"))
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 7.12.1. Plotting the effect sizes ----
# ----------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_dp.dead_estimates$parameter_plot,
                     season = lions_output_dp.dead_estimates$season,
                     mean = lions_output_dp.dead_estimates$mean,
                     median = lions_output_dp.dead_estimates$median,
                     BCI90_upper = lions_output_dp.dead_estimates$upperCI,
                     BCI90_lower = lions_output_dp.dead_estimates$lowerCI,
                     overlapping0 = FALSE)

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean detection\nprobability"))
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "DP_Dead.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

dpdead_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior, # Density plot
                       aes(x = posterior, y = variable, fill = season), 
                       color = NA, size = 1, alpha = 0.3, 
                       position = position_dodge(width = 0.7)) +
  geom_linerange(data = df.plot, # Credible intervals
                 aes(xmin = BCI90_lower, xmax = BCI90_upper, y = variable, 
                     color = season), 
                 position = position_dodge(width = 0.7), size = 1, linetype = "solid", show.legend = F) +
  geom_point(data = df.plot, # Median
             aes(x = median, y = variable, color = season), 
             position = position_dodge(width = 0.7), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", # 0-line
             color = colPlot, size = 0.4) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-10, 10, 1), limits = c(-3.9, 1)) +
  theme_general() +
  ggtitle("Death recovery probability") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

dpdead_plot

dev.off()


## 7.12.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
dp.dead.estimate = function(season = "wet"){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.dp = lions_output_GLMM[, grep("mu.dp.dead", 
                                       colnames(lions_output_GLMM))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_dp.dead.epsilons[, seq(1, 59, 2)], 1, median)

  }
  
  else{
    
    # Get dry-season parameter values
    mean.dp = lions_output_GLMM[, grep("mu.dp.dead", 
                                       colnames(lions_output_GLMM))][, 2]
    
    # Get dry-season epsilons
    median.epsilons= apply(lions_output_dp.dead.epsilons[, seq(2, 60, 2)], 1, median)

  }
  
  # Calculate vital-rate prediction
  pred =  mean.dp + median.epsilons
  
  return(pred)
  
}


# Average detection probability per season

# Empty prediction dataframe
dp.dead.season.estimates = expand.grid(season = c("wet", "dry"))

# Fill in prediction and credible intervals
dp.dead.season.estimates$pred = plogis(apply(apply(dp.dead.season.estimates,
                                                   1, 
                                                   FUN = function(x) dp.dead.estimate(season = x[1])), 
                                             2, median))
dp.dead.season.estimates$lwr = plogis(apply(apply(dp.dead.season.estimates, 
                                                  1, 
                                                  FUN = function(x) dp.dead.estimate(season = x[1])), 
                                            2, FUN = function(x) quantile(x, probs = 0.05)))
dp.dead.season.estimates$upr = plogis(apply(apply(dp.dead.season.estimates,
                                                  1, 
                                                  FUN = function(x) dp.dead.estimate(season = x[1])),
                                            2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_DP_Dead_Season.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(dp.dead.season.estimates, aes(x = season, y = pred)) +
  geom_point(position = position_dodge(0.8), colour = colPlot) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5, colour = colPlot) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  xlab("Season") +
  ylab("Death recovery probability") +
  labs(title = "Death recovery probability") +
  theme_minimal() +
  theme(axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        panel.grid = element_blank(),
        axis.text = element_text(colour = colPlot, size = 8, family = font),
        axis.title = element_text(colour = colPlot, size = 9, family = font),
        legend.title = element_text(colour = colPlot, size = 9, family = font),
        legend.text = element_text(colour = colPlot, size = 8, family = font),
        plot.title = element_blank(),
        legend.position = "none")

dev.off()




###########################################################################
#
# 8. Epsilon posterior distributions ----
#
###########################################################################

# Dataframe to match full parameter names, model output column, and file name
# We add two empty slots for plot aesthetics reasons.

params.labels = data.frame(# Output model column name 
                           col.label = c("s.sa1",
                                         "s.sa2",
                                         "s.af",
                                         "s.ym",
                                         "s.nm",
                                         "s.rm",
                                         "emig.ym",
                                         "t.ym.nm",
                                         "takeover",
                                         "eviction",
                                         "dp.pride",
                                         "dp.nm",
                                         "dp.dead",
                                         "1",
                                         "2"),
                           # Full parameter name
                           param.name = c("Young-subadult survival",
                                          "Old-subadult survival",
                                          "Adult-female survival",
                                          "Young-male survival",
                                          "Nomadic-male survival",
                                          "Resident-male survival",
                                          "Young-male emigration",
                                          "Young- to nomadic-male\ntransition",
                                          "Nomadic-male takeover",
                                          "Resident-male eviction",
                                          "Pride detection probability",
                                          "Nomadic-male\ndetection probability",
                                          "Dead recovery",
                                          "1",
                                          "2"),
                           # File name
                           param.file = c("SA1Surv",
                                          "SA2Surv",
                                          "AFSurv",
                                          "YMSurv",
                                          "NMSurv",
                                          "RMSurv",
                                          "YMEmig",
                                          "YMNMTrans",
                                          "Takeover",
                                          "Eviction",
                                          "PrideDP",
                                          "NMDP",
                                          "DeadDP",
                                          "1",
                                          "2"))

# Epsilon distribution dataframe
nSamples = 1

epsilon_df = data.frame(parameter = rep(rep(params.labels$param.name, 
                                            each = nSamples * length(seq(1985, 2014))), 2),
                        label = rep(rep(params.labels$col.label,
                                        each = nSamples * length(seq(1985, 2014))), 2),
                        season = rep(c("Wet", "Dry"), 
                                     each = nSamples * length(seq(1985, 2014)) * length(params.labels$param.name)),
                        year = rep(rep(rep(seq(1985, 2014), each = nSamples), 
                                       length(params.labels$param.name)), 2), 
                        median = NA, 
                        lwr = NA,
                        upr = NA)


# Fill in epsilon distribution dataframe
for(p in 1:nrow(params.labels[1:13, ])){
  
  param = params.labels[p, ] # Get parameter full name
  
  epsilon.colname = paste0("epsilon.", param[1]) # Get model output column name
  
  # Subset wet-season espilons
  epsilons_wet = lions_output_GLMM[, grep(epsilon.colname, 
                                          colnames(lions_output_GLMM))] 
  epsilons_wet = epsilons_wet[, seq(1, 59, 2)]
  
  # Get mean epsilons and credible intervals
  epsilon_df$median[which(epsilon_df$label == as.character(param[1]) &
                          epsilon_df$season == "Wet")] = apply(epsilons_wet, 2,
                                                               median)
  epsilon_df$upr[which(epsilon_df$label == as.character(param[1]) &
                       epsilon_df$season == "Wet")] = apply(epsilons_wet, 2, 
                                                            FUN = function(x) quantile(x, probs = 0.95))
  epsilon_df$lwr[which(epsilon_df$label == as.character(param[1]) &
                       epsilon_df$season == "Wet")] = apply(epsilons_wet, 2, 
                                                              FUN = function(x) quantile(x, probs = 0.05))
  
  # Subset dry-season espilons
  epsilons_dry = lions_output_GLMM[, grep(epsilon.colname, colnames(lions_output_GLMM))] 
  epsilons_dry = epsilons_dry[, seq(2, 60, 2)]
  
  # Get mean epsilons and credible intervals
  epsilon_df$median[which(epsilon_df$label == as.character(param[1]) &
                          epsilon_df$season == "Dry")] = apply(epsilons_dry, 2, 
                                                               median)
  epsilon_df$upr[which(epsilon_df$label == as.character(param[1]) &
                       epsilon_df$season == "Dry")] = apply(epsilons_dry, 2, 
                                                              FUN = function(x) quantile(x, probs = 0.95))
  epsilon_df$lwr[which(epsilon_df$label == as.character(param[1]) &
                       epsilon_df$season == "Dry")] = apply(epsilons_dry, 2, 
                                                              FUN = function(x) quantile(x, probs = 0.05))
}

# Factorize parameter names
epsilon_df$parameter = factor(epsilon_df$parameter,
                              levels = unique(epsilon_df$parameter))


# Plot yearly epsilons
plot.epsilon.post = 
  ggplot(epsilon_df, aes(x = year, y = median, colour = season)) +
  facet_wrap(~ parameter, ncol = 3, scales = "free_y") +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) + 
  scale_fill_manual(name = "Season",
                    values = c(cbbPalette[4], cbbPalette[2]),
                    labels = c("Wet", "Dry")) +
  scale_colour_manual(name = "Season",
                      values = c(cbbPalette[4], cbbPalette[2]),
                      labels = c("Wet", "Dry")) +
  scale_y_continuous(n.breaks = 5) +
  xlab("Year") +
  ylab("Epsilon (year random effect)") + 
  theme_minimal() %+replace%    
  theme(axis.text = element_text(colour = colPlot, size = 7, family = font),
        axis.title.x = element_text(colour = colPlot, size = 8, family = font, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = 8, family = font, margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, linewidth = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, family = font, hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "right", 
        legend.justification = "right",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(colour = colPlot, size = 7, family = font),
        legend.title = element_text(colour = colPlot, size = 8, family = font),
        legend.key.size = unit(15, "pt"))

png(filename = paste0("EpsilonPosteriorDistributions.png"), 
    width = 16, 
    height = 17, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

print(plot.epsilon.post)

dev.off()