############################################################################
#
# This script uses samples obtained from the chains of an MCMC algorithm.
#
# The aim of this script is to simulate new capture-history datasets 
# based on sampled parameters to perform posterior predictive checks.
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

library(coda)


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
# 3. Simulating datasets ----
#
###########################################################################

## 3.1. Sampling posterior distributions ----
# --------------------------------------

# 500 samples from the posterior distribution for each parameter
posterior_sampled = lions_output_GLMM[seq(1, nrow(lions_output_GLMM), 
                                          length.out = 500), ]


## 3.2. Simulating datasets ----
# -------------------------

# Empty simulated capture histories
simulated_ch = lions.ch
simulated_ch[, ] = NA

# Add state at first capture
for(i in 1:nrow(simulated_ch)){
  
  simulated_ch[i, lions.first[i]] = lions.ch[i, lions.first[i]] 
  
}


# Function to fill in the simulated capture histories
# using the observed covariate values and predicted vital rates
transition_function = function(# Initial values and datasets
                               initial_capture, # First capture occasion of a lion 
                               stages_true,     # True stages
                               stages_observed, # Observed stages
                               
                               
                               # Intercepts of vital-rate models
                               mu.s.sa1,        # Young-subadult survival
                               mu.s.sa2f,       # Female old-subadult survival
                               mu.s.sa2m,       # Male old-subadult survival
                               mu.s.af,         # Adult-female survival
                               mu.s.ym,         # Young-male survival
                               mu.s.nm,         # Nomadic-male survival
                               mu.s.rm,         # Resident-male survival
                               mu.emig.ym,      # Young-male emigration 
                               mu.t.ym.nm,      # Young-male transition to nomadic male
                               mu.takeover,     # Nomadic-male takeover
                               mu.eviction,     # Resident-male eviction
                               mu.dp.pride,     # Pride-individual detection
                               mu.dp.nm,        # Nomadic-male detection
                               mu.dp.dead,      # Dead-individual detection
                               
                               # Betas of vital-rate models
                               # Young-subadult survival
                               s.sa1.beta.nb.nm.coal.hr,       # Number of nomadic coalitions in the home range
                               s.sa1.beta.nb.af.pride,         # Number of adult females in the pride
                               s.sa1.beta.habitat.woodland,    # Habitat (NA for grassland which is the reference level)
                               # Old-subadult survival
                               s.sa2.beta.nb.nm.coal.hr,       # Number of nomadic coalitions in the home range 
                               s.sa2.beta.nb.af.pride,         # Number of adult females in the pride
                               s.sa2.beta.habitat.woodland,    # Habitat (NA for grassland which is the reference level)
                               # Adult-female survival
                               s.af.beta.age,                  # Age
                               s.af.beta.nb.af.pride,          # Number of adult females in the pride
                               s.af.beta.nb.nm.coal.hr,        # Number of nomadic coalitions in the home range
                               s.af.beta.habitat.woodland,     # Habitat (NA for grassland which is the reference level)
                               # Young-male survival
                               s.ym.beta.nb.nm.coal.hr,        # Number of nomadic coalitions in the home range
                               s.ym.beta.nb.af.pride,          # Number of adult females in the pride
                               s.ym.beta.habitat.woodland,     # Habitat (NA for grassland which is the reference level)
                               # Nomadic-male survival
                               s.nm.beta.coal.size,            # Coalition size
                               s.nm.beta.habitat.woodland,     # Habitat (NA for grassland which is the reference level)
                               # Resident-male survival
                               s.rm.beta.coal.size,            # Coalition size
                               s.rm.beta.nb.nm.coal.hr,        # Number of nomadic coalitions in the home range
                               s.rm.beta.habitat.woodland,     # Habitat (NA for grassland which is the reference level)
                               # Young-male emigration
                               emig.ym.beta.habitat.woodland,  # Habitat (NA for grassland which is the reference level) 
                               # Nomadic-male takeover
                               takeover.beta.coal.size,        # Coalition size
                               takeover.beta.habitat.woodland, # Habitat (NA for grassland which is the reference level) 
                               # Resident-male eviction
                               eviction.beta.coal.size,        # Coalition size
                               eviction.beta.nb.nm.coal.hr,    # Number of nomadic coalitions in the home range
                               eviction.beta.habitat.woodland, # Habitat (NA for grassland which is the reference level) 
                               # Pride-individual detection
                               dp.pride.beta.habitat.woodland, # Habitat (NA for grassland which is the reference level) 
                               # Nomadic-male detection
                               dp.nm.beta.habitat.woodland,    # Habitat (NA for grassland which is the reference level) 
                               
                               # Random effect epsilons
                               epsilon.s.sa1,     # Young-subadult survival  
                               epsilon.s.sa2,     # Old-subadult survival
                               epsilon.s.af,      # Adult-female survival
                               epsilon.s.ym,      # Young-male survival
                               epsilon.s.nm,      # Nomadic-male survival
                               epsilon.s.rm,      # Resident-male survival
                               epsilon.emig.ym,   # Young-male emigration
                               epsilon.t.ym.nm,   # Young-male transition to nomadic male
                               epsilon.takeover,  # Nomadic-male takeover
                               epsilon.eviction,  # Resident-male eviction
                               epsilon.dp.pride,  # Pride-individual detection
                               epsilon.dp.nm,     # Nomadic-male detection
                               epsilon.dp.dead,   # Dead-individual detection
                               
                               # Covariates
                               groups,          # Pride or male coalition to which a lion belongs
                               seasons,         # Season at a given timestep
                               years,           # Year at a given timestep
                               nb.nm.coal.hr_matrix, # Unscaled number of nomadic coalitions in the home range
                               # Parameters for covariate sampling and standardization
                               min.nb.nm.coal.hr = min(nb.nm.coal.hr.unscaled, na.rm = T), # Minimum observed number of nomadic coalitions in the home range
                               max.nb.nm.coal.hr = max(nb.nm.coal.hr.unscaled, na.rm = T), # Maximum observed number of nomadic coalitions in the home range
                               mu.nb.nm.coal.hr = mean(nb.nm.coal.hr.unscaled, na.rm = T), # Mean observed number of nomadic coalitions in the home range
                               sd.nb.nm.coal.hr = sd(nb.nm.coal.hr.unscaled, na.rm = T),   # Observed standard deviation of number of nomadic coalitions in the home range
                               nb.nm.coal.hr.theta = nb.nm.coal.hr.theta,                  # Theta of the negative binomial model for the number of nomadic coalitions in the home range
                               nb.af.pride_matrix, # Unscaled number of adult females in the pride
                               # Parameters for covariate sampling and standardization
                               min.nb.af.pride = min(nb.af.pride.unscaled, na.rm = T), # Minimum observed number of adult females in the pride
                               max.nb.af.pride = max(nb.af.pride.unscaled, na.rm = T), # Maximum observed number of adult females in the pride
                               mu.nb.af.pride = mean(nb.af.pride.unscaled, na.rm = T), # Observed mean number of adult females in the pride
                               sd.nb.af.pride = sd(nb.af.pride.unscaled, na.rm = T),   # Observed standard deviation of the number of adult females in the pride
                               nb.af.pride.theta = nb.af.pride.theta,                  # Theta of the negative binomial model for the nnumber of adult females in the pride
                               habitat_vector, # Observed habitat covariate
                               habitat.prob = habitat.prob, # Probability to be in the plains
                               age_vector, # Individual age at each timestep
                               coal.size_matrix, # Unscaled coalition size
                               # Parameters for covariate sampling and standardization
                               min.coal.size = min(coal.size.unscaled, na.rm = T), # Minimum observed coalition size
                               max.coal.size = max(coal.size.unscaled, na.rm = T), # Maximum observed coalition size
                               mu.coal.size = mean(coal.size.unscaled, na.rm = T), # Observed mean coalition size
                               sd.coal.size = sd(coal.size.unscaled, na.rm = T),   # Observed standard deviation of coalition size
                               coal.size.lambda = coal.size.lambda){               # Lambda of the Poisson model for coalition size
  
  
  # Initial stage and capture
  stage_true = stages_true[initial_capture]
  stage_observed = stages_observed[initial_capture]
  
  # For each capture occasion
  for(capture in initial_capture:(length(stages_true) - 1)){
    
    # Sampling covariate values and standardizing
    # We use truncated distributions to sample the missing values to stay
    # within the observed values.
    # Number of nomadic coalitions in the home range
    nb.nm.coal.hr = ifelse(!is.na(nb.nm.coal.hr_matrix[groups[capture], capture]), 
                           nb.nm.coal.hr_matrix[groups[capture], capture], 
                           (truncdist::rtrunc(1, "nbinom", 
                                              a = min.nb.nm.coal.hr, 
                                              b = max.nb.nm.coal.hr, 
                                              prob = 0.35, 
                                              size = nb.nm.coal.hr.theta) 
                            - mu.nb.nm.coal.hr) / (2 * sd.nb.nm.coal.hr))
    
    # Number of adult females in the pride
    nb.af.pride = ifelse(!is.na(nb.af.pride_matrix[groups[capture], capture]), 
                         nb.af.pride_matrix[groups[capture], capture], 
                         (truncdist::rtrunc(1, "nbinom", 
                                            a = min.nb.af.pride,
                                            b = max.nb.af.pride, 
                                            prob = 0.35, 
                                            size = nb.af.pride.theta)
                          - mu.nb.af.pride) / (2 * sd.nb.af.pride))
    
    # Habitat
    habitat = ifelse(!is.na(habitat_vector[capture]), 
                     habitat_vector[capture], 
                     rbinom(1, size = 1, prob = habitat.prob) + 1)
    
    # Age
    age = age_vector[capture]
    
    # Coalition size
    coal.size = ifelse(!is.na(coal.size_matrix[groups[capture], capture]), 
                       coal.size_matrix[groups[capture], capture], 
                       (truncdist::rtrunc(1, "pois",
                                          a = min.coal.size,
                                          b = max.coal.size, 
                                          lambda = coal.size.lambda)
                        - mu.coal.size) / (2 * sd.coal.size))
    
    
    # Calculate vital-rate predictions
    
    # Young-subadult survival
    survSA1 = as.numeric(plogis(mu.s.sa1[seasons[capture]] + 
                                         s.sa1.beta.nb.nm.coal.hr[seasons[capture]] * nb.nm.coal.hr +
                                         s.sa1.beta.nb.af.pride[seasons[capture]] * nb.af.pride +
                                         s.sa1.beta.habitat.woodland[seasons[capture], habitat] +
                                         epsilon.s.sa1[seasons[capture], years[capture]]))
    
    # Female old-subadult survival
    survSA2F = as.numeric(plogis(mu.s.sa2f[seasons[capture]] + 
                                          s.sa2.beta.nb.nm.coal.hr * nb.nm.coal.hr +
                                          s.sa2.beta.nb.af.pride[seasons[capture]] * nb.af.pride +
                                          s.sa2.beta.habitat.woodland[seasons[capture], habitat] +
                                          epsilon.s.sa2[seasons[capture], years[capture]]))
                                  
    # Male old-subadult survival       
    survSA2M = as.numeric(plogis(mu.s.sa2m[seasons[capture]] + 
                                          s.sa2.beta.nb.nm.coal.hr * nb.nm.coal.hr +
                                          s.sa2.beta.nb.af.pride[seasons[capture]] * nb.af.pride +
                                          s.sa2.beta.habitat.woodland[seasons[capture], habitat] +
                                          epsilon.s.sa2[seasons[capture], years[capture]]))
    
    # Adult-female survival  
    survAF = as.numeric(plogis(mu.s.af[seasons[capture]] + 
                                        s.af.beta.age[seasons[capture]] * age +
                                        s.af.beta.nb.nm.coal.hr[seasons[capture]] * nb.nm.coal.hr +
                                        s.af.beta.nb.af.pride[seasons[capture]] * nb.af.pride +
                                        s.af.beta.habitat.woodland[seasons[capture], habitat] +
                                        epsilon.s.af[seasons[capture], years[capture]]))
    
    # Young-male survival
    survYM = as.numeric(plogis(mu.s.ym[seasons[capture]] + 
                                        s.ym.beta.nb.nm.coal.hr[seasons[capture]] * nb.nm.coal.hr +
                                        s.ym.beta.nb.af.pride[seasons[capture]] * nb.af.pride +
                                        s.ym.beta.habitat.woodland[seasons[capture], habitat] +
                                        epsilon.s.ym[seasons[capture], years[capture]]))
    
    # Nomadic-male survival
    survNM = as.numeric(plogis(mu.s.nm[seasons[capture]] + 
                                        s.nm.beta.coal.size[seasons[capture]] * coal.size +
                                        s.nm.beta.habitat.woodland[seasons[capture], habitat] +
                                        epsilon.s.nm[seasons[capture], years[capture]]))
    
    # Resident-male survival
    survRM = as.numeric(plogis(mu.s.rm[seasons[capture]] + 
                                        s.rm.beta.nb.nm.coal.hr[seasons[capture]] * nb.nm.coal.hr +
                                        s.rm.beta.coal.size[seasons[capture]] * coal.size +
                                        s.rm.beta.habitat.woodland[seasons[capture], habitat] +
                                        epsilon.s.rm[seasons[capture], years[capture]]))
    
    # Young-male emigration
    emigYM = as.numeric(plogis(mu.emig.ym[seasons[capture]] + 
                                        emig.ym.beta.habitat.woodland[seasons[capture], habitat] +
                                        epsilon.emig.ym[seasons[capture], years[capture]]))
    
    # Young-male transition to nomadic male
    transYM = as.numeric(plogis(mu.t.ym.nm[seasons[capture]] + 
                                         epsilon.t.ym.nm[seasons[capture], years[capture]]))
    
    # Nomadic-male takeover
    takeover = as.numeric(plogis(mu.takeover[seasons[capture]] + 
                                          takeover.beta.coal.size[seasons[capture]] * coal.size +
                                          takeover.beta.habitat.woodland[seasons[capture], habitat] +
                                          epsilon.takeover[seasons[capture], years[capture]]))
    
    # Resident-male eviction
    eviction = as.numeric(plogis(mu.eviction[seasons[capture]] + 
                                          eviction.beta.coal.size[seasons[capture]] * coal.size +
                                          eviction.beta.nb.nm.coal.hr[seasons[capture]] * nb.nm.coal.hr +
                                          eviction.beta.habitat.woodland[seasons[capture], habitat] +
                                          epsilon.eviction[seasons[capture], years[capture]]))
    
    # Pride-individual detection
    dpPride = as.numeric(plogis(mu.dp.pride[seasons[capture + 1]] +
                                         dp.pride.beta.habitat.woodland[seasons[capture + 1], habitat] +
                                         epsilon.dp.pride[seasons[capture + 1], years[capture + 1]]))
    
    # Nomadic-male detection
    dpNM = as.numeric(plogis(mu.dp.nm[seasons[capture + 1]] +
                                      dp.nm.beta.habitat.woodland[seasons[capture + 1], habitat] +
                                      epsilon.dp.nm[seasons[capture + 1], years[capture + 1]]))
    
    # Dead-individual detection
    dpDead = as.numeric(plogis(mu.dp.dead[seasons[capture + 1]] +
                                        epsilon.dp.dead[seasons[capture + 1], years[capture + 1]]))
    
    
    # Transitions between true stages
    # SA1 = Young subadult
    # SA2F/SA2M = Female and male old subadult
    # AF = Adult female
    # YM = Young male
    # NM = Nomadic male
    # RM = Resident male
    
    # SA1-SA2F/SA2M
    if(stage_true == 1){
    
      stage_true = rbinom(1, 1, prob = survSA1) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 2 * rbinom(1, 1, prob = 0.55) # Does the individual go to SA2F?
        if(stage_true == 0) stage_true = 3 # If not, the individual goes to SA2M.
      }
    }
    
    # SA2F-AF
    else if(stage_true == 2){stage_true = 4 * rbinom(1, 1, prob = survSA2F)} # Does the individual survive?
    
    # SA2M-YM1
    else if(stage_true == 3){stage_true = 5 * rbinom(1, 1, prob = survSA2M)} # Does the individual survive?
    
    # AF-AF
    else if(stage_true == 4){stage_true = 4 * rbinom(1, 1, prob = survAF)}   # Does the individual survive?
    
    # YM1-YM2/NM/RM
    else if(stage_true == 5){
      
      stage_true = rbinom(1, 1, prob = survYM) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 6 * rbinom(1, 1, prob = (1 - emigYM)) # Does the individual stay in its pride?
        
        if(stage_true == 0){ # The individual leaves its pride.
          
          stage_true = 9 * rbinom(1, 1, prob = transYM) # Does the individual become nomadic?
          if(stage_true == 0) stage_true = 10 # If not, the individual becomes resident.
        }
      }
    }
    
    # YM2-YM3/NM/RM
    else if(stage_true == 6){
      
      stage_true = rbinom(1, 1, prob = survYM) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 7 * rbinom(1, 1, prob = (1 - emigYM)) # Does the individual stay in its pride?
        
        if(stage_true == 0){ # The individual leaves its pride.
          
          stage_true = 9 * rbinom(1, 1, prob = transYM) # Does the individual become nomadic?
          if(stage_true == 0) stage_true = 10 # If not, the individual becomes resident.
        }
      }
    } 
    
    # YM3-YM4/NM/RM
    else if(stage_true == 7){
      
      stage_true = rbinom(1, 1, prob = survYM) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 8 * rbinom(1, 1, prob = (1 - emigYM)) # Does the individual stay in its pride?
        
        if(stage_true == 0){ # The individual leaves its pride.
          
          stage_true = 9 * rbinom(1, 1, prob = transYM) # Does the individual become nomadic?
          if(stage_true == 0) stage_true = 10 # If not, the individual becomes resident.
        }
      }
    }
    
    # YM4-NM/RM
    else if(stage_true == 8){
      
      stage_true = rbinom(1, 1, prob = survYM) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 9 * rbinom(1, 1, prob = transYM) # Does the individual become nomadic?
        
        if(stage_true == 0) stage_true = 10 # If not, the individual becomes resident.
      }
    }
    
    # NM-NM/RM
    else if(stage_true == 9){
      
      stage_true = rbinom(1, 1, prob = survNM) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 9 * rbinom(1, 1, prob = (1 - takeover)) # Does the individual stays nomadic?
        
        if(stage_true == 0) stage_true = 10 # If not, the individual becomes resident.
      }
    }
    
    # RM-RM/NM
    else if(stage_true == 10){
      
      stage_true = rbinom(1, 1, prob = survRM) # Does the individual survive?
      
      if(stage_true != 0){ # The individual survives.
        
        stage_true = 10 * rbinom(1, 1, prob = (1 - eviction)) # Does the individual stays resident?
        
        if(stage_true == 0) stage_true = 9 # If not, the individual becomes nomadic.
      }
    }
    
    # Add current true stage to true stage data
    stages_true[capture + 1] = stage_true
    
    
    # Assess observation using detection probabilities
    
    # Pride individual
    if(stage_true %in% c(1:8, 10)){
      
      stage_observed = stage_true * rbinom(1, 1, prob = dpPride) # Is the individual observed ?
      
      # If yes, updated stage observed with true stage, otherwise assign unobserved stage
      if(stage_observed == stage_true) stages_observed[capture + 1] = stage_observed
      else stages_observed[capture + 1] = 13
      
    }
    
    # Nomadic male
    else if(stage_true == 9){
      
      stage_observed = stage_true * rbinom(1, 1, prob = dpNM) # Is the individual observed ?
      
      # If yes, updated stage observed with true stage, otherwise assign unobserved stage
      if(stage_observed == stage_true) stages_observed[capture + 1] = stage_observed
      else stages_observed[capture + 1] = 13
      
    }
    
    # Alive-Newly dead-Permanently dead
    if(stage_true == 0){ # If the individual dies
      
      # Next capture stage is dead, the followings are permanently dead
      stage_true = 11
      stages_true[capture + 1] = stage_true
      stages_true[(capture + 2):length(stages_true)] = 12
      
      
      # Assess observation of dead individual
      stage_observed = stage_true * rbinom(1, 1, prob = dpDead) # Is the individual observed ?
      
      # If yes, updated stage observed with true stage, otherwise assign unobserved stage
      if(stage_observed == stage_true) stages_observed[capture + 1] = stage_observed
      else stages_observed[capture + 1] = 13
      
      stages_observed[capture + 2:length(stages_observed)] = 13
      
      break
    }
  }
  
  return(stages_observed)
}



# Simulate datasets in parallel

ncpus = 5 # Number of cores

posterior_sampled_slices = trunc(seq(0, nrow(posterior_sampled), # Posterior samples for each core
                                     length.out = ncpus + 1))


# Function to give to the parallelization function to simulate the data
simulate_data = function(sim){ # sim = core
  
  if (sim > 1) {      # Give each core a random seed
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
  }
  
  # Get posterior samples for the core
  posterior_sampled_parallel = posterior_sampled[(posterior_sampled_slices[sim] + 1):posterior_sampled_slices[sim + 1], ]
  
  # List to store simulated capture history for each set of sampled parameters
  list_simulated_ch = vector(mode = "list",
                             length = 10 * nrow(posterior_sampled_parallel))
  
  k = 1 # Initialize counter to store data in the list
   
  # For each set of sampled posterior parameter distribution
  for(post in 1:nrow(posterior_sampled_parallel)){
    
    print(paste("Posterior samples set: ", post))
    
    # Simulate 10 datasets per sample
    for(sim_dat in 1:10){
      
      print(paste("Dataset: ", sim_dat))
      
      simulated_ch = lions.ch[0, ] # Empty capture histories
      
      # For each lion
      for(i in 1:nrow(lions.ch)){
        
        # Simulate capture history 
        simulated_ch = rbind(simulated_ch,
        transition_function(# Initial values and datasets
                            initial_capture = as.numeric(lions.first[i]), # First capture occasion of a lion 
                            stages_true = lions.ch[i, ],                  # True stages
                            stages_observed = lions.ch[i, ],              # Observed stages
                            
                            # Intercepts of vital-rate models
                            mu.s.sa1 = posterior_sampled_parallel[post,    # Young-subadult survival
                                                                  grep("mu.s.sa1", 
                                                                       colnames(posterior_sampled_parallel))],
                            mu.s.sa2f = posterior_sampled_parallel[post,   # Female old-subadult survival
                                                                   grep("mu.s.sa2f", 
                                                                        colnames(posterior_sampled_parallel))],
                            mu.s.sa2m = posterior_sampled_parallel[post,   # Male old-subadult survival
                                                                   grep("mu.s.sa2m", 
                                                                        colnames(posterior_sampled_parallel))],
                            mu.s.af = posterior_sampled_parallel[post,     # Adult-female survival
                                                                 grep("mu.s.af", 
                                                                      colnames(posterior_sampled_parallel))],
                            mu.s.ym = posterior_sampled_parallel[post,     # Young-male survival
                                                                 grep("mu.s.ym", 
                                                                      colnames(posterior_sampled_parallel))],
                            mu.s.nm = posterior_sampled_parallel[post,     # Nomadic-male survival
                                                                 grep("mu.s.nm", 
                                                                      colnames(posterior_sampled_parallel))],
                            mu.s.rm = posterior_sampled_parallel[post,     # Resident-male survival
                                                                 grep("mu.s.rm", 
                                                                      colnames(posterior_sampled_parallel))],
                            mu.emig.ym = posterior_sampled_parallel[post,  # Young-male emigration
                                                                    grep("mu.emig.ym", 
                                                                         colnames(posterior_sampled_parallel))],
                            mu.t.ym.nm = posterior_sampled_parallel[post,  # Young-male transition to nomadic male
                                                                    grep("mu.t.ym.nm",
                                                                         colnames(posterior_sampled_parallel))],
                            mu.takeover = posterior_sampled_parallel[post, # Nomadic-male takeover
                                                                     grep("mu.takeover", 
                                                                          colnames(posterior_sampled_parallel))],
                            mu.eviction = posterior_sampled_parallel[post, # Resident-male eviction
                                                                     grep("mu.eviction", 
                                                                          colnames(posterior_sampled_parallel))],
                            mu.dp.pride = posterior_sampled_parallel[post, # Pride-individual detection
                                                                     grep("mu.dp.pride", 
                                                                          colnames(posterior_sampled_parallel))],
                            mu.dp.nm = posterior_sampled_parallel[post,    # Nomadic-male detection
                                                                  grep("mu.dp.nm", 
                                                                       colnames(posterior_sampled_parallel))],
                            mu.dp.dead = posterior_sampled_parallel[post,  # Dead-individual detection
                                                                    grep("mu.dp.dead",
                                                                         colnames(posterior_sampled_parallel))],
                            
                            # Betas of vital-rate models
                            # Young-subadult survival
                            s.sa1.beta.nb.nm.coal.hr = posterior_sampled_parallel[post,              # Number of nomadic coalitions in the home range
                                                                                  grep("s.sa1.beta.nb.nm.coal.hr", 
                                                                                       colnames(posterior_sampled_parallel))],
                            s.sa1.beta.nb.af.pride = posterior_sampled_parallel[post,                # Number of adult females in the pride
                                                                                grep("s.sa1.beta.nb.af.pride", 
                                                                                     colnames(posterior_sampled_parallel))],
                            s.sa1.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,    # Habitat (NA for grassland which is the reference level)
                                                                                            grep("s.sa1.beta.habitat.woodland", 
                                                                                                 colnames(posterior_sampled_parallel))], 2, 2),
                            # Old-subadult survival
                            s.sa2.beta.nb.nm.coal.hr = posterior_sampled_parallel[post,              # Number of nomadic coalitions in the home range 
                                                                                  grep("s.sa2.beta.nb.nm.coal.hr", 
                                                                                       colnames(posterior_sampled_parallel))],
                            s.sa2.beta.nb.af.pride = posterior_sampled_parallel[post,                # Number of adult females in the pride
                                                                                grep("s.sa2.beta.nb.af.pride", 
                                                                                     colnames(posterior_sampled_parallel))],
                            s.sa2.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,    # Habitat (NA for grassland which is the reference level)
                                                                                            grep("s.sa2.beta.habitat.woodland",
                                                                                                 colnames(posterior_sampled_parallel))], 2, 2),
                            # Adult-female survival
                            s.af.beta.age = posterior_sampled_parallel[post,                         # Age
                                                                       grep("s.af.beta.age",
                                                                            colnames(posterior_sampled_parallel))],
                            s.af.beta.nb.af.pride = posterior_sampled_parallel[post,                 # Number of adult females in the pride
                                                                               grep("s.af.beta.nb.af.pride", 
                                                                                    colnames(posterior_sampled_parallel))],
                            s.af.beta.nb.nm.coal.hr = posterior_sampled_parallel[post,               # Number of nomadic coalitions in the home range
                                                                                 grep("s.af.beta.nb.nm.coal.hr",
                                                                                      colnames(posterior_sampled_parallel))],
                            s.af.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,     # Habitat (NA for grassland which is the reference level)
                                                                                           grep("s.af.beta.habitat.woodland", 
                                                                                                colnames(posterior_sampled_parallel))], 2, 2),
                            # Young-male survival
                            s.ym.beta.nb.nm.coal.hr = posterior_sampled_parallel[post,               # Number of nomadic coalitions in the home range
                                                                                 grep("s.ym.beta.nb.nm.coal.hr", 
                                                                                      colnames(posterior_sampled_parallel))],
                            s.ym.beta.nb.af.pride = posterior_sampled_parallel[post,                 # Number of adult females in the pride
                                                                               grep("s.ym.beta.nb.af.pride", 
                                                                                    colnames(posterior_sampled_parallel))],
                            s.ym.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,     # Habitat (NA for grassland which is the reference level)
                                                                                           grep("s.ym.beta.habitat.woodland",
                                                                                                colnames(posterior_sampled_parallel))], 2, 2),
                            # Nomadic-male survival
                            s.nm.beta.coal.size = posterior_sampled_parallel[post,                   # Coalition size
                                                                             grep("s.nm.beta.coal.size",
                                                                                  colnames(posterior_sampled_parallel))],
                            s.nm.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,     # Habitat (NA for grassland which is the reference level)
                                                                                           grep("s.nm.beta.habitat.woodland",
                                                                                                colnames(posterior_sampled_parallel))], 2, 2),
                            # Resident-male survival
                            s.rm.beta.coal.size = posterior_sampled_parallel[post,                   # Coalition size
                                                                             grep("s.rm.beta.coal.size", 
                                                                                  colnames(posterior_sampled_parallel))],
                            s.rm.beta.nb.nm.coal.hr = posterior_sampled_parallel[post,               # Number of nomadic coalitions in the home range
                                                                                 grep("s.rm.beta.nb.nm.coal.hr", 
                                                                                      colnames(posterior_sampled_parallel))],
                            s.rm.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,     # Habitat (NA for grassland which is the reference level)
                                                                                           grep("s.rm.beta.habitat.woodland", 
                                                                                                colnames(posterior_sampled_parallel))], 2, 2),
                            # Young-male emigration
                            emig.ym.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,  # Habitat (NA for grassland which is the reference level) 
                                                                                              grep("emig.ym.beta.habitat.woodland", 
                                                                                                   colnames(posterior_sampled_parallel))], 2, 2),
                            # Nomadic-male takeover
                            takeover.beta.coal.size = posterior_sampled_parallel[post,               # Coalition size
                                                                                 grep("takeover.beta.coal.size", 
                                                                                      colnames(posterior_sampled_parallel))],
                            takeover.beta.habitat.woodland = matrix(posterior_sampled_parallel[post, # Habitat (NA for grassland which is the reference level) 
                                                                                               grep("takeover.beta.habitat.woodland", 
                                                                                                    colnames(posterior_sampled_parallel))], 2, 2),
                            # Resident-male eviction
                            eviction.beta.coal.size = posterior_sampled_parallel[post,               # Coalition size
                                                                                 grep("eviction.beta.coal.size", 
                                                                                      colnames(posterior_sampled_parallel))],
                            eviction.beta.nb.nm.coal.hr = posterior_sampled_parallel[post,           # Number of nomadic coalitions in the home range
                                                                                     grep("eviction.beta.nb.nm.coal.hr", 
                                                                                          colnames(posterior_sampled_parallel))],
                            eviction.beta.habitat.woodland = matrix(posterior_sampled_parallel[post, # Habitat (NA for grassland which is the reference level) 
                                                                                               grep("eviction.beta.habitat.woodland",
                                                                                                    colnames(posterior_sampled_parallel))], 2, 2),
                            # Pride-individual detection
                            dp.pride.beta.habitat.woodland = matrix(posterior_sampled_parallel[post, # Habitat (NA for grassland which is the reference level) 
                                                                                               grep("dp.pride.beta.habitat.woodland",
                                                                                                    colnames(posterior_sampled_parallel))], 2, 2),
                            # Nomadic-male detection
                            dp.nm.beta.habitat.woodland = matrix(posterior_sampled_parallel[post,    # Habitat (NA for grassland which is the reference level) 
                                                                                            grep("dp.nm.beta.habitat.woodland", 
                                                                                                 colnames(posterior_sampled_parallel))], 2, 2),
                            
                            # Random effect epsilons
                            epsilon.s.sa1 = matrix(posterior_sampled_parallel[post,    # Young-subadult survival
                                                                              grep("epsilon.s.sa1",
                                                                                   colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.s.sa2 = matrix(posterior_sampled_parallel[post,    # Old-subadult survival
                                                                              grep("epsilon.s.sa2", 
                                                                                   colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.s.af = matrix(posterior_sampled_parallel[post,     # Adult-female survival
                                                                             grep("epsilon.s.af", 
                                                                                  colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.s.ym = matrix(posterior_sampled_parallel[post,     # Young-male survival
                                                                             grep("epsilon.s.ym",
                                                                                  colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.s.nm = matrix(posterior_sampled_parallel[post,     # Nomadic-male survival
                                                                             grep("epsilon.s.nm", 
                                                                                  colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.s.rm = matrix(posterior_sampled_parallel[post,     # Resident-male survival
                                                                             grep("epsilon.s.rm", 
                                                                                  colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.emig.ym = matrix(posterior_sampled_parallel[post,  # Young-male emigration
                                                                                grep("epsilon.emig.ym",
                                                                                     colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.t.ym.nm = matrix(posterior_sampled_parallel[post,  # Young-male transition to nomadic male
                                                                                grep("epsilon.t.ym.nm",
                                                                                     colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.takeover = matrix(posterior_sampled_parallel[post, # Nomadic-male takeover
                                                                                 grep("epsilon.takeover", 
                                                                                      colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.eviction = matrix(posterior_sampled_parallel[post, # Resident-male eviction
                                                                                 grep("epsilon.eviction",
                                                                                      colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.dp.pride = matrix(posterior_sampled_parallel[post, # Pride-individual detection
                                                                                 grep("epsilon.dp.pride", 
                                                                                      colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.dp.nm = matrix(posterior_sampled_parallel[post,    # Nomadic-male detection
                                                                              grep("epsilon.dp.nm",
                                                                                   colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            epsilon.dp.dead = matrix(posterior_sampled_parallel[post,  # Dead-individual detection
                                                                                grep("epsilon.dp.dead",
                                                                                     colnames(posterior_sampled_parallel))], 2, length(unique(year))),
                            
                            # Covariates
                            groups = lions.groups[i, ], # Pride or male coalition to which a lion belongs
                            seasons = season,           # Season at a given timestep
                            years = year,               # Year at a given timestep
                            nb.nm.coal.hr_matrix = nb.nm.coal.hr, # Unscaled number of nomadic coalitions in the home range
                            # Parameters for covariate sampling and standardization 
                            min.nb.nm.coal.hr = min(nb.nm.coal.hr.unscaled, na.rm = T), # Minimum observed number of nomadic coalitions in the home range
                            max.nb.nm.coal.hr = max(nb.nm.coal.hr.unscaled, na.rm = T), # Maximum observed number of nomadic coalitions in the home range
                            mu.nb.nm.coal.hr = mean(nb.nm.coal.hr.unscaled, na.rm = T), # Mean observed number of nomadic coalitions in the home range
                            sd.nb.nm.coal.hr = sd(nb.nm.coal.hr.unscaled, na.rm = T),   # Observed standard deviation of number of nomadic coalitions in the home range
                            nb.nm.coal.hr.theta = nb.nm.coal.hr.theta,                  # Theta of the negative binomial model for the number of nomadic coalitions in the home range
                            nb.af.pride_matrix = nb.af.pride, # Unscaled number of adult females in the pride
                            # Parameters for covariate sampling and standardization
                            min.nb.af.pride = min(nb.af.pride.unscaled, na.rm = T),     # Minimum observed number of adult females in the pride
                            max.nb.af.pride = max(nb.af.pride.unscaled, na.rm = T),     # Maximum observed number of adult females in the pride
                            mu.nb.af.pride = mean(nb.af.pride.unscaled, na.rm = T),     # Observed mean number of adult females in the pride
                            sd.nb.af.pride = sd(nb.af.pride.unscaled, na.rm = T),       # Observed standard deviation of the number of adult females in the pride
                            nb.af.pride.theta = nb.af.pride.theta,                      # Theta of the negative binomial model for the nnumber of adult females in the pride
                            habitat_vector = habitat[i, ], # Observed habitat covariate
                            habitat.prob = habitat.prob,   # Probability to be in the plains
                            age_vector = age[i, ],         # Individual age at each timestep
                            coal.size_matrix = coal.size,  # Unscaled coalition size
                            # Parameters for covariate sampling and standardization
                            min.coal.size = min(coal.size.unscaled, na.rm = T), # Minimum observed coalition size
                            max.coal.size = max(coal.size.unscaled, na.rm = T), # Maximum observed coalition size
                            mu.coal.size = mean(coal.size.unscaled, na.rm = T), # Observed mean coalition size
                            sd.coal.size = sd(coal.size.unscaled, na.rm = T),   # Observed standard deviation of coalition size
                            coal.size.lambda = coal.size.lambda))               # Lambda of the Poisson model for coalition size

      }
      
      list_simulated_ch[[k]] = simulated_ch # Add simulated capture histories to the list
      
      k = k + 1
      
    }
  }
  
  return(list_simulated_ch)
}


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run.

# Set up parallel environment
library(snowfall)

sfInit(parallel = TRUE, cpus = ncpus, # Initialisation and progress output file
       slaveOutfile = "DataSimulationProgress.txt")
sfExport(list = c(ls(), ".Random.seed")) # Export current environment to each core
sfLibrary("snowfall", character.only = TRUE) # Load snowfall 

# Results of each set of iterations are saved to a file
data_simulation_output = sfClusterApplyLB(1:ncpus, simulate_data) # Run simulations on each core

save(data_simulation_output, file = "Simulated_CH.RData")
sfStop() # Stop parallelization