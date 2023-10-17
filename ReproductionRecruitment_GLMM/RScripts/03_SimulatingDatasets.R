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
lions.data = read.csv("Data/01_LionsDemographicData.csv")

# MCMC samples
load("LionFullModel_Repro_Rec_Output.RData")

# Full output as a matrix
lions_output_fullGLMM = as.matrix(rbind(lions_results_fullGLMM[[1]],
                                        lions_results_fullGLMM[[2]],
                                        lions_results_fullGLMM[[3]],
                                        lions_results_fullGLMM[[4]])) 




###########################################################################
#
# 2. Formatting dataset ----
#
###########################################################################

# Subset female data only and remove nomadic females
females.data = lions.data[which(lions.data$stage == "AF"), ]
females.data = females.data[- which(females.data$pride == "NO"), ]


# Format age and habitat
lions.data$age.at.capture = lions.data$age.at.capture / 12 # Age in years
lions.data$habitat.code = lions.data$habitat.code + 1 # Habitat as c(1, 2) instead of c(0, 1)

# Get probability of being in the woodland for missing covariate values
habitat = females.data$habitat.code - 1
summary(glm(c(habitat) ~ 1, "binomial"))
habitat.intercept.estimate = coef(glm(c(habitat) ~ 1, "binomial"))
barplot(table(habitat)/sum(table(habitat)), ylim = c(0, 1))
points(dbinom(0:1, size = 1, prob = plogis(habitat.intercept.estimate)))
habitat.prob = as.numeric(plogis(habitat.intercept.estimate))


# New season number column
females.data$season.nb = NA
females.data$season.nb[which(females.data$season == "wet")] = 1
females.data$season.nb[which(females.data$season == "dry")] = 2


# New year number column
year.number = data.frame(cbind(unique(females.data$year), # Assigning a number to each year of the dataset
                               seq(1:length(unique(females.data$year))))) 
colnames(year.number) = c("year","year.nb")


females.data$year.nb = NA

for(row in 1:nrow(females.data)){ # Add number corresponding to each year in the data
  
  females.data$year.nb[row] = year.number$year.nb[which(year.number$year == females.data$year[row])]
  
}


# Standardize covariates

females.data$age.at.capture.scaled = (females.data$age.at.capture - mean(females.data$age.at.capture, na.rm = T)) / 
                                     (2 * sd(females.data$age.at.capture, na.rm = T))
females.data$nb.af.pride.scaled = (females.data$nb_af_pride - mean(females.data$nb_af_pride, na.rm = T)) /
                                     (2 * sd(females.data$nb_af_pride, na.rm = T))
females.data$nb.nm.coal.hr.scaled = (females.data$nb_nm_coal_hr - mean(females.data$nb_nm_coal_hr, na.rm = T)) / 
                                    (2 * sd(females.data$nb_nm_coal_hr, na.rm = T))






###########################################################################
#
# 3. Simulating datasets ----
#
###########################################################################


## 3.1. Sampling posterior distributions ----
# --------------------------------------

# 500 samples from the posterior distribution for each parameter
posterior_sampled = lions_output_fullGLMM[seq(1, nrow(lions_output_fullGLMM), 
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
list_simulated_data = vector(mode = "list", 
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
    
    list_simulated_data_seasonal_coeffs[[k]] = simulated_data # Add simulated data to the list
    
    k = k + 1
    
  }
}


# Save simulated data list
save(list_simulated_data, file = "Repro_Rec_Simulated_Data.RData")