############################################################################
#
# This script uses the capture-recapture data of a lion population in the 
# Serengeti between 1984 and 2014.
# The aim of this script is to model the reproduction probability and 
# recruitment of lions in the population using two GLMMs 
# (binomial and Poisson).
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

library(nimble)
library(coda)
library(ggpubr)
library(boot)
library(bayestestR)

nimbleOptions(disallow_multivariate_argument_expressions = F)


## 1.3. Loading data ----
# ------------------

# Lion demographic dataset
females.data = read.csv("Data/01_LionsFemalesDemographicData.csv")




###########################################################################
#
# 2. Format data ----
#
###########################################################################

# Standardize covariates
females.data$age.at.capture.scaled = (females.data$age.at.capture - mean(females.data$age.at.capture, na.rm = T)) /
                                     (2 * sd(females.data$age.at.capture, na.rm = T))
females.data$nb.af.pride.scaled = (females.data$nb_af_pride - mean(females.data$nb_af_pride, na.rm = T)) / 
                                  (2 * sd(females.data$nb_af_pride, na.rm = T))
females.data$nb.nm.coal.hr.scaled = (females.data$nb_nm_coal_hr - mean(females.data$nb_nm_coal_hr, na.rm = T)) / 
                                    (2 * sd(females.data$nb_nm_coal_hr, na.rm = T))


# Get habitat distribution
summary(glm((females.data$habitat.code) - 1 ~ 1, "binomial"))
habitat.intercept.estimate = coef(glm((females.data$habitat.code) - 1 ~ 1, "binomial"))
barplot(table((females.data$habitat.code) - 1) / sum(table((females.data$habitat.code) - 1)), 
        ylim = c(0, 1))
points(dbinom(0:1, size = 1, prob = inv.logit(habitat.intercept.estimate))) # Probability of being in the woodland




###########################################################################
#
# 3. Preparing and fitting the model ----
#
###########################################################################

## 3.1. Preparing the model ----
# -------------------------

# Model code 

lions_code_recruit_repro = nimbleCode({
  # -------------------------------------------------
  #
  # Parameters:
  # 
  # mu.repro = mean reproduction probability 
  #
  # mu.rec = mean recruitment 
  # -------------------------------------------------
  #
  # Linear models
  
  for(row in 1:R){ # For each row of the dataset
    
    # Reproduction
    
    reproduction[row] ~ dbern(p[row]) # Reproduction is a binary response 
                                      # sampled from a Bernoulli distribution
                                      # with a mean p that depends on the estimated parameters
    logit(p[row]) <- mu.repro[season[row]] + 
                     repro.beta.habitat.woodland[season[row], habitat[row]] +
                     repro.beta.age[season[row]] * age[row] +
                     repro.beta.quad.age[season[row]] * (age[row] * age[row]) +
                     repro.beta.nb.af.pride[season[row]] * nb.af.pride[row] +
                     repro.beta.quad.nb.af.pride[season[row]] * (nb.af.pride[row] * nb.af.pride[row]) +
                     repro.beta.nb.nm.coal.hr[season[row]] * nb.nm.coal.hr[row] + 
                     repro.beta.nb.af.pride.nb.nm.coal.hr[season[row]] * nb.af.pride[row] * nb.nm.coal.hr[row] +
                     epsilon.repro[season[row], year[row]]
    
    
    # Recruitment
    
    cubs[row] ~ dpois(lambda[row]) # Recruitment is an integer
                                   # sampled from a Poisson distribution
                                   # with a mean lambda that depends on the estimated parameters
    log(lambda[row]) <- mu.rec[season[row]] + 
                        rec.beta.habitat.woodland[season[row], habitat[row]] +
                        rec.beta.nb.af.pride[season[row]] * nb.af.pride[row] +
                        rec.beta.nb.nm.coal.hr[season[row]] * nb.nm.coal.hr[row] + 
                        rec.beta.nb.af.pride.nb.nm.coal.hr[season[row]] * nb.af.pride[row] * nb.nm.coal.hr[row] +
                        epsilon.rec[season[row], year[row]]
    
  } 
  
  
  # Priors and constraints
  
  # Fixed effects 
  
  for(seas in 1:2){
    
    # Means
    
    mu.repro[seas] ~ dunif(-10, 10)                             # Reproduction probability           
    mu.rec[seas] ~ dunif(-10, 2)                                # Recruitment
    
    
    # Other parameters
    
    # Reproduction probability
    
    repro.beta.habitat.woodland[seas, 1] <- 0                   # Habitat - grassland (reference level)
    repro.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)       # Habitat - woodland
    
    repro.beta.age[seas] ~ dunif(-10, 10)                       # Age
    repro.beta.quad.age[seas] ~ dunif(-10, 10)                  # Quadratic age
    repro.beta.nb.af.pride[seas] ~ dunif(-10, 10)               # Number of adult females in the pride
    repro.beta.quad.nb.af.pride[seas] ~ dunif(-10, 10)          # Quadratic number of adult females in the pride
    repro.beta.nb.nm.coal.hr[seas] ~ dunif(-10, 10)             # Number of nomadic coalitions in the home range
    repro.beta.nb.af.pride.nb.nm.coal.hr[seas] ~ dunif(-10, 10) # Interaction between the number of
                                                                # adult females in the pride and the
                                                                # number of nomadic coalitions in 
                                                                # the home range
    
    
    # Recruitment
    
    rec.beta.habitat.woodland[seas, 1] <- 0                     # Habitat - grassland (reference level)
    rec.beta.habitat.woodland[seas, 2] ~ dunif(-5, 1)           # Habitat - woodland
    
    rec.beta.nb.af.pride[seas] ~ dunif(-5, 1)                   # Number of adult females in the pride
    rec.beta.nb.nm.coal.hr[seas] ~ dunif(-5, 1)                 # Number of nomadic coalitions in the home range
    rec.beta.nb.af.pride.nb.nm.coal.hr[seas] ~ dunif(-5, 1)     # Interaction between the number of
                                                                # adult females in the pride and the
                                                                # number of nomadic coalitions in 
                                                                # the home range
  }
  
  
  # Season-specific yearly random effects
  
  for(seas in 1:2){
    
    sigma.repro[seas] ~ dunif(0, 10) # Reproduction probability
    sigma.rec[seas] ~ dunif(0, 10)   # Recruitment
    
    for(yr in 1:30){
      
      epsilon.repro[seas, yr] ~ dnorm(0, sd = sigma.repro[seas]) # Reproduction probability
      epsilon.rec[seas, yr] ~ dnorm(0, sd = sigma.rec[seas])     # Recruitment
      
    }
    
  }
  
  # Sampling of missing habitat covariates values using the observed distribution
  
  for(row in 1:R){
    
    habitat[row] <- temp.habitat[row] + 1
    temp.habitat[row] ~ dbin(size = 1, prob = habitat.prob)
  }
})


# Model constants

lions_constants_GLMM <- list(R = nrow(females.data), # Number of rows in the data
                             season = females.data$season.nb, # Season
                             year = females.data$year.nb,     # Year
                             age = females.data$age.at.capture.scaled, # Age
                             nb.af.pride = females.data$nb.af.pride.scaled,     # Number of adult females in the pride
                             nb.nm.coal.hr = females.data$nb.nm.coal.hr.scaled, # Number of nomadic coalitions in the home range 
                             habitat.prob = inv.logit(habitat.intercept.estimate)) # Observed probability of being in the woodland habitat


# Model data

lions_data_fullGLMM <- list(reproduction = females.data$reproduction,     # Reproduction 
                            cubs = females.data$cubs,                     # Recruitment
                            habitat = females.data$habitat.code,          # Habitat as a (1, 2) variable for the model
                            temp.habitat = females.data$habitat.code - 1) # Habitat as a (0, 1) variable for the sampling of missing covariate values


# Model initial values

# Initial values for the sampling of missing habitat values
inits.habitat = function(cov){ # cov is the observed data
  
  cov.new = cov
  which.na = which(is.na(cov))
  
  # Replace NAs by sampling the observed negative binomial distribution
  if(length(which.na) != 0){
    
    cov.new[which.na] = rbinom(length(which.na), size = 1, prob = 0.3)
    cov.new[-which.na] = NA
  }
  
  return(cov.new) # Full covariate data
}


# Setting initial values for the model parameters and creating the list
# of initial values

lions_inits_fullGLMM <- function() {
  list(# Intercepts
       mu.repro = runif(2, -2, 2), # Reproduction probability
       mu.rec = runif(2, -1, 1),   # Recruitment
       
       # Betas
       # Reproduction
       repro.beta.habitat.woodland = matrix(c(rep(NA, 2),            # Habitat (NA for grassland which is the reference level)
                                              runif(2, 
                                                    -0.05, 
                                                    0.05)), 
                                            2, 2),
       repro.beta.age = runif(2, -0.05, 0.05),                       # Age
       repro.beta.quad.age = runif(2, -0.05, 0.05),                  # Quadratic age
       repro.beta.nb.af.pride = runif(2, -0.05, 0.05),               # Number of adult females in the pride
       repro.beta.quad.nb.af.pride = runif(2, -0.05, 0.05),          # Quadratic number of adult females in the pride
       repro.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),             # Number of nomadic coalitions in the home range
       repro.beta.nb.af.pride.nb.nm.coal.hr = runif(2, -0.05, 0.05), # Interaction between the number of
                                                                     # adult females in the pride and the
                                                                     # number of nomadic coalitions in 
                                                                     # the home range 
       # Recruitment
       rec.beta.habitat.woodland = matrix(c(rep(NA, 2),              # Habitat (NA for grassland which is the reference level)
                                            runif(2,
                                                  -0.05, 
                                                  0.05)), 
                                          2, 2),
       rec.beta.nb.af.pride = runif(2, -0.05, 0.05),                 # Number of adult females in the pride
       rec.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),               # Number of nomadic coalitions in the home range
       rec.beta.nb.af.pride.nb.nm.coal.hr = runif(2, -0.05, 0.05),   # Interaction between the number of
                                                                     # adult females in the pride and the
                                                                     # number of nomadic coalitions in 
                                                                     # the home range 
       
       # Random effect sigmas
       sigma.repro = runif(2, 0, 10), # Reproduction probability
       sigma.rec = runif(2, 0, 10),   # Recruitment
       
       # Random effect epsilons
       epsilon.repro = matrix(0,      # Reproduction probability
                              nrow = length(unique(females.data$season.nb)), 
                              ncol = length(unique(females.data$year.nb))),
       epsilon.rec = matrix(0,        # Recruitment
                            nrow = length(unique(females.data$season.nb)),
                            ncol = length(unique(females.data$year.nb))),
       
       # Missing habitat initial values
       temp.habitat = inits.habitat(females.data$habitat.code - 1))
}


# Sample initial values
init.val = list(lions_inits_fullGLMM(), # Chain 1
                lions_inits_fullGLMM(), # Chain 2
                lions_inits_fullGLMM(), # Chain 3
                lions_inits_fullGLMM()) # Chain 4


# Parameters monitored

parameters = c(# Intercepts
               "mu.repro",                             # Reproduction probability
               "mu.rec",                               # Recruitment
               
               # Betas
               # Reproduction probability
               "repro.beta.habitat.woodland",          # Habitat
               "repro.beta.age",                       # Age
               "repro.beta.quad.age",                  # Quadratic age
               "repro.beta.nb.af.pride",               # Number of adult females in the pride
               "repro.beta.quad.nb.af.pride",          # Quadratic number of adult females in the pride
               "repro.beta.nb.nm.coal.hr",             # Number of nomadic coalitions in the home range
               "repro.beta.nb.af.pride.nb.nm.coal.hr", # Interaction between the number of
                                                       # adult females in the pride and the
                                                       # number of nomadic coalitions in 
                                                       # the home range 
               # Recruitment 
               "rec.beta.habitat.woodland",            # Habitat
               "rec.beta.nb.af.pride",                 # Number of adult females in the pride
               "rec.beta.nb.nm.coal.hr",               # Number of nomadic coalitions in the home range
               "rec.beta.nb.af.pride.nb.nm.coal.hr",   # Interaction between the number of
                                                       # adult females in the pride and the
                                                       # number of nomadic coalitions in 
                                                       # the home range 
               
               # Random effect sigmas
               "sigma.repro",                          # Reproduction probability
               "sigma.rec",                            # Recruitment
               
               # Random effect epsilons
               "epsilon.repro",                        # Reproduction probability
               "epsilon.rec")                          # Recruitment


## 3.2. Fitting the model ----
# -----------------------

lions_results_fullGLMM = 
  nimbleMCMC(code = lions_code_recruit_repro,
             constants = lions_constants_GLMM,
             data = lions_data_fullGLMM,
             inits = lions_inits_fullGLMM,
             monitors = parameters,
             niter = 55000, 
             nburnin = 10000,
             nchains = 4,
             thin = 1,
             samplesAsCodaMCMC = TRUE)

save(lions_results_fullGLMM, file = "FullModel_Repro_Rec.RData")
