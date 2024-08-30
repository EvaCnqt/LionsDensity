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

rm(list = ls(all.names = T))

## 1.1. Loading libraries ----
# -----------------------


## 1.2. Loading data ----
# ------------------

# Lion demographic dataset
females.data = read.csv("Data/01_LionsFemalesDemographicData.csv")

# MCMC samples
load("Output/Lions_Reproduction_Recruitment_MCMCSamples.RData")


# Density dependent covariates for covariate standardization

# Number of females in a pride
nb.af.pride.unscaled = read.csv("Data/023_Covariate_NbAFpride.csv", 
                                stringsAsFactors = F, row.names = 1)
nb.af.pride.unscaled = as.matrix(nb.af.pride.unscaled)
range(nb.af.pride.unscaled, na.rm = T)

min.nb.af.pride = min(nb.af.pride.unscaled, na.rm = T)
max.nb.af.pride = max(nb.af.pride.unscaled, na.rm = T)
mu.nb.af.pride = mean(nb.af.pride.unscaled, na.rm = T)
sd.nb.af.pride = sd(nb.af.pride.unscaled, na.rm = T)

# Age
age.unscaled = read.csv("Data/024_Covariate_Age.csv", row.names = 1)
age.unscaled = as.matrix(age.unscaled)
range(age.unscaled, na.rm = T)
mu.age = mean(age.unscaled, na.rm = T)
sd.age = sd(age.unscaled, na.rm = T)

# Number of nomadic coalitions in the home range of a pride or 
# a resident male coalition
nb.nm.coal.hr.unscaled = read.csv("Data/025_Covariate_NbNMCoalHR.csv", row.names = 1)
nb.nm.coal.hr.unscaled = as.matrix(nb.nm.coal.hr.unscaled)
range(nb.nm.coal.hr.unscaled, na.rm = T)
min.nb.nm.coal.hr = min(nb.nm.coal.hr.unscaled, na.rm = T)
max.nb.nm.coal.hr = max(nb.nm.coal.hr.unscaled, na.rm = T)
mu.nb.nm.coal.hr = mean(nb.nm.coal.hr.unscaled, na.rm = T)
sd.nb.nm.coal.hr = sd(nb.nm.coal.hr.unscaled, na.rm = T)




###########################################################################
#
# 2. Format data ----
#
###########################################################################

# Standardize covariates
females.data$age.at.capture.scaled = (females.data$age.at.capture - mu.age) /
  (2 * sd.age)
females.data$nb.af.pride.scaled = (females.data$nb_af_pride - mu.nb.af.pride) / 
  (2 * sd.nb.af.pride)
females.data$nb.nm.coal.hr.scaled = (females.data$nb_nm_coal_hr - mu.nb.nm.coal.hr) / 
  (2 * sd.nb.nm.coal.hr)


# Get habitat proportion
habitat = females.data$habitat.code - 1
summary(glm(c(habitat) ~ 1, "binomial"))
habitat.intercept.estimate = coef(glm(c(habitat) ~ 1, "binomial"))
barplot(table(habitat)/sum(table(habitat)), ylim = c(0, 1))
points(dbinom(0:1, size = 1, prob = boot::inv.logit(habitat.intercept.estimate)))
habitat.prob = as.numeric(boot::inv.logit(habitat.intercept.estimate))




###########################################################################
#
# 3. Simulating datasets ----
#
###########################################################################


## 3.1. Sampling posterior distributions ----
# --------------------------------------

# 500 samples from the posterior distribution for each parameter
posterior_sampled = lions_output_repro_recruit[seq(1, nrow(lions_output_repro_recruit), 
                                              length.out = 500), ]


## 3.2. Simulating datasets ----
# -------------------------

# Function to create the simulated data
# using the observed covariate values and predicted vital rates
repro_rec_function = function(data = females.data, # Observed data
                              
                              # Intercepts of vital-rate models
                              mu.repro, # Reproduction probability
                              mu.rec,   # Recruitment
                              
                              # Betas of vital-rate models
                              # Reproduction probability
                              repro.beta.habitat.woodland,          # Habitat
                              repro.beta.age,                       # Age
                              repro.beta.quad.age,                  # Quadratic age
                              repro.beta.nb.af.pride,               # Number of adult females in the pride
                              repro.beta.quad.nb.af.pride,          # Quadratic number of adult females in the pride
                              repro.beta.nb.nm.coal.hr,             # Number of nomadic coalitions in the home range
                              repro.beta.nb.af.pride.nb.nm.coal.hr, # Interaction between the number of
                                                                    # adult females in the pride and the
                                                                    # number of nomadic coalitions in 
                                                                    # the home range
                              # Recruitment
                              rec.beta.habitat.woodland,            # Habitat
                              rec.beta.nb.af.pride,                 # Number of adult females in the pride
                              rec.beta.nb.nm.coal.hr,               # Number of nomadic coalitions in the home range
                              rec.beta.nb.af.pride.nb.nm.coal.hr,   # Interaction between the number of
                                                                    # adult females in the pride and the
                                                                    # number of nomadic coalitions in 
                                                                    # the home range
                              
                              # Random effect epsilons
                              epsilon.repro,  # Reproduction probability
                              epsilon.rec,    # Recruitment
                              
                              # Covariates
                              seasons = females.data$season.nb, # Season at a given timestep
                              years = females.data$year.nb,     # Year at a given timestep
                              nb.nm.coal.hr_vector = females.data$nb.nm.coal.hr.scaled, # Number of nomadic coalitions in the home range
                              nb.af.pride_vector = females.data$nb.af.pride.scaled, # Number of adult females in the pride
                              habitat_vector = females.data$habitat.code, # Habitat
                              habitat.prob = habitat.prob, # Habitat probability for missing habitat values
                              age_vector = females.data$age.at.capture.scaled){ # Age
  
  
  # Remove values of reproduction and recruitment
  data$reproduction = NA
  data$cubs = NA
  
  # For each data row
  for(i in 1:nrow(data)){
    
    # Getting covariate values and sampling missing habitat values 
    habitat = ifelse(!is.na(habitat_vector[i]), habitat_vector[i], rbinom(1, size = 1, prob = habitat.prob) + 1)
    
    age = age_vector[i]
    
    nb.nm.coal.hr = nb.nm.coal.hr_vector[i]
    
    nb.af.pride = nb.af.pride_vector[i]
    
    
    # Calculate vital-rate predictions
    
    # Reproduction probability
    repro.prob = as.numeric(plogis(mu.repro[seasons[i]] + 
                                   repro.beta.habitat.woodland[seasons[i], habitat] +
                                   repro.beta.age[seasons[i]] * age +
                                   repro.beta.quad.age[seasons[i]] * age^2 +
                                   repro.beta.nb.af.pride[seasons[i]] * nb.af.pride +
                                   repro.beta.quad.nb.af.pride[seasons[i]] * nb.af.pride^2 +
                                   repro.beta.nb.nm.coal.hr[seasons[i]] * nb.nm.coal.hr +
                                   repro.beta.nb.af.pride.nb.nm.coal.hr[seasons[i]] * nb.af.pride * nb.nm.coal.hr +
                                   epsilon.repro[seasons[i], years[i]]))
    
    # Recruitment
    cubs.mean = as.numeric(exp(mu.rec[seasons[i]] + 
                               rec.beta.habitat.woodland[seasons[i], habitat] +
                               rec.beta.nb.af.pride[seasons[i]] * nb.af.pride +
                               rec.beta.nb.nm.coal.hr[seasons[i]] * nb.nm.coal.hr +
                               rec.beta.nb.af.pride.nb.nm.coal.hr[seasons[i]] * nb.af.pride * nb.nm.coal.hr +
                               epsilon.rec[seasons[i], years[i]]))
    
    
    # Predicted reproduction and recruitment values
    data$reproduction[i] = rbinom(1, 1, prob = repro.prob)
    data$cubs[i] = rpois(1, lambda = cubs.mean)
    
  }
  
  return(data)
}


# Empty list to store the simulated datasets
lions_repro_recruit_simulated_data = vector(mode = "list", 
                                            length = 10 * nrow(posterior_sampled))

k = 1 # Initialize counter to store data in the list


# For each set of sampled posterior parameter distribution
for(post in 1:nrow(posterior_sampled)){
  
  print(paste("Posterior samples set: ", post))
  
  # Simulate 10 datasets per sample
  for(sim_dat in 1:10){
    
    print(paste("Dataset: ", sim_dat))
    
    # Simulate data
    simulated_data = repro_rec_function(data = females.data, # Observed data
                                        
                                        # Intercepts of vital-rate models
                                        mu.repro = posterior_sampled[post,   # Reproduction probability
                                                                     grep("mu.repro",
                                                                          colnames(posterior_sampled))],
                                        mu.rec = posterior_sampled[post,     # Recruitment
                                                                   grep("mu.rec", 
                                                                        colnames(posterior_sampled))],
                                        
                                        # Betas of vital-rate models
                                        # Reproduction probability
                                        repro.beta.habitat.woodland = matrix(posterior_sampled[post,   # Habitat
                                                                                               grep("repro.beta.habitat.woodland", 
                                                                                                    colnames(posterior_sampled))],
                                                                             2, 2),
                                        repro.beta.age = posterior_sampled[post,                       # Age
                                                                           grep("repro.beta.age", 
                                                                                colnames(posterior_sampled))],
                                        repro.beta.quad.age = posterior_sampled[post,                  # Quadratic age
                                                                                grep("repro.beta.quad.age", 
                                                                                     colnames(posterior_sampled))],
                                        repro.beta.nb.af.pride = posterior_sampled[post,               # Number of adult females in the pride
                                                                                   grep("repro.beta.nb.af.pride", 
                                                                                        colnames(posterior_sampled))][c(1, 2)],
                                        repro.beta.quad.nb.af.pride = posterior_sampled[post,          # Quadratic number of adult females in the pride
                                                                                        grep("repro.beta.quad.nb.af.pride",
                                                                                             colnames(posterior_sampled))],
                                        repro.beta.nb.nm.coal.hr = posterior_sampled[post,             # Number of nomadic coalitions in the home range
                                                                                     grep("repro.beta.nb.nm.coal.hr", 
                                                                                          colnames(posterior_sampled))],
                                        repro.beta.nb.af.pride.nb.nm.coal.hr = posterior_sampled[post, # Interaction between the number of
                                                                                                       # adult females in the pride and the
                                                                                                       # number of nomadic coalitions in 
                                                                                                       # the home range
                                                                                                 grep("repro.beta.nb.af.pride.nb.nm.coal.hr",
                                                                                                      colnames(posterior_sampled))],
                                        # Recruitment
                                        rec.beta.habitat.woodland = matrix(posterior_sampled[post,     # Habitat
                                                                                             grep("rec.beta.habitat.woodland", 
                                                                                                  colnames(posterior_sampled))], 
                                                                           2, 2),
                                        rec.beta.nb.af.pride = posterior_sampled[post,                 # Number of adult females in the pride
                                                                                 grep("rec.beta.nb.af.pride",
                                                                                      colnames(posterior_sampled))][c(1, 2)],
                                        rec.beta.nb.nm.coal.hr = posterior_sampled[post,               # Number of nomadic coalitions in the home range
                                                                                   grep("rec.beta.nb.nm.coal.hr", 
                                                                                        colnames(posterior_sampled))],
                                        rec.beta.nb.af.pride.nb.nm.coal.hr = posterior_sampled[post,   # Interaction between the number of
                                                                                                       # adult females in the pride and the
                                                                                                       # number of nomadic coalitions in 
                                                                                                       # the home range
                                                                                               grep("rec.beta.nb.af.pride.nb.nm.coal.hr", 
                                                                                                    colnames(posterior_sampled))],
                                        
                                        # Random effect epsilons
                                        epsilon.repro = matrix(posterior_sampled[post, # Reproduction probability
                                                                                 grep("epsilon.repro",
                                                                                      colnames(posterior_sampled))],
                                                               2, length(unique(females.data$year))),
                                        epsilon.rec = matrix(posterior_sampled[post,   # Recruitment
                                                                               grep("epsilon.rec",
                                                                                    colnames(posterior_sampled))], 
                                                             2, length(unique(females.data$year))),
                                        
                                        # Covariates
                                        seasons = females.data$season.nb, # Season at a given timestep
                                        years = females.data$year.nb,     # Year at a given timestep
                                        nb.nm.coal.hr_vector = females.data$nb.nm.coal.hr.scaled, # Number of nomadic coalitions in the home range
                                        nb.af.pride_vector = females.data$nb.af.pride.scaled,     # Number of adult females in the pride
                                        habitat_vector = females.data$habitat.code, # Habitat
                                        habitat.prob = habitat.prob, # Habitat probability for missing habitat values
                                        age_vector = females.data$age.at.capture.scaled) # Age
    
    lions_repro_recruit_simulated_data[[k]] = simulated_data # Add simulated data to the list
    
    k = k + 1
    
  }
}


# Save simulated data list
save(lions_repro_recruit_simulated_data, file = "Output/Lions_Reproduction_Recruitment_Simulated_Data.RData")
