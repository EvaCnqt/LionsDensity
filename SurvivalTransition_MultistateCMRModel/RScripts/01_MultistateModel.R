############################################################################
#
# The aim of this script is to fit a multistate model to estimate the
# stage-specific survival and transition rates, as well as detection 
# probabilities in a population of African lions from the Serengeti National
# Park, Tanzania. The data used was collected between 1984 and 2014.
# This model also aims at understanding the effect of environmental (habitat
# and rainfall) and density-dependent factors (number of females in the pride, 
# male coalition size, or number of nomadic coalitions) on vital rates.
#
# Author: Eva Conquet
#
############################################################################


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
library(boot)

nimbleOptions(disallow_multivariate_argument_expressions = F)


## 1.3. Loading data ----
# ------------------

# Custom distribution
source("RScripts/00_dDHMM_lionKF.R")

# Individual capture histories
lions.ch = read.csv("Data/01_LionsCaptureHistories.csv", row.names = 1)
lions.ch = as.matrix(lions.ch)


## 1.4. Covariates and their observed distribution ----
# ------------------------------------------------

# The observed distribution will be used to sample the missing values
# of covariates needed in the model.


# Year and season
year = read.csv("Data/021_Covariate_Year.csv", stringsAsFactors = F)$x
year = rep(c(1:30), each = 2)

season = read.csv("Data/022_Covariate_Season.csv", stringsAsFactors = F)$x


# Habitat
habitat = read.csv("Data/023_Covariate_Habitat.csv", row.names = 1)
habitat = as.matrix(habitat)

# Distribution
summary(glm(c(habitat) ~ 1, "binomial"))
habitat.intercept.estimate = coef(glm(c(habitat) ~ 1, "binomial"))
barplot(table(habitat) / sum(table(habitat)), ylim = c(0, 1))
points(dbinom(0:1, size = 1, prob = inv.logit(habitat.intercept.estimate)))
habitat = habitat + 1 # To avoid zeros in the NIMBLE model


# Density-dependent covariates

# Number of females in a pride
nb.af.pride.unscaled = read.csv("Data/024_Covariate_NbAFpride.csv", 
                                stringsAsFactors = F, 
                                row.names = 1)
nb.af.pride.unscaled = as.matrix(nb.af.pride.unscaled)
range(nb.af.pride.unscaled, na.rm = T)
nb.af.pride = (nb.af.pride.unscaled - mean(nb.af.pride.unscaled, na.rm = T)) / 
  (2 * sd(nb.af.pride.unscaled, na.rm = T)) # Standardize covariate

# Distribution
summary(MASS::glm.nb(c(nb.af.pride.unscaled) ~ 1))
negbin.model = MASS::glm.nb(c(nb.af.pride.unscaled) ~ 1)
hist(nb.af.pride.unscaled, freq = F)
lines(dnbinom(min(nb.af.pride.unscaled, na.rm = T):max(nb.af.pride.unscaled, na.rm = T), 
              9.7, 
              mu = mean(nb.af.pride.unscaled, na.rm = T)))
hist(rnbinom(10000, 
             size = negbin_model$theta,
             mu = mean(nb.af.pride.unscaled, na.rm = T)), 
     freq = F, col = "red", add = T)
hist(rnbinom(10000, 
             size = negbin_model$theta, prob = 0.65), 
     freq = F, col = "red", add = T)
nb.af.pride.theta = negbin.model$theta


# Age
age.unscaled = read.csv("Data/025_Covariate_Age.csv", row.names = 1)
age.unscaled = as.matrix(age.unscaled)
range(age.unscaled, na.rm = T)
age = (age.unscaled - mean(age.unscaled, na.rm = T)) / 
  (2 * sd(age.unscaled, na.rm = T)) # Standardize covariate


# Male coalition size 
coal.size.unscaled = read.csv("Data/026_Covariate_CoalSize.csv", row.names = 1)
coal.size.unscaled = as.matrix(coal.size.unscaled)
range(coal.size.unscaled, na.rm = T)
coal.size = (coal.size.unscaled - mean(coal.size.unscaled, na.rm = T)) / 
  (2 * sd(coal.size.unscaled, na.rm = T)) # Standardize covariate

# Distribution
summary(glm(c(coal.size.unscaled) ~ 1, "poisson"))
summary(glm(c(coal.size.unscaled) ~ 1, "quasipoisson"))
poisson.model = glm(c(coal.size.unscaled) ~ 1, "quasipoisson")
summary(MASS::glm.nb(c(coal.size.unscaled) ~ 1))

hist(coal.size.unscaled, freq = F)

# Poisson
lines(dpois(min(coal.size.unscaled, na.rm = T):max(coal.size.unscaled, na.rm = T), 2.35))
hist(rpois(10000, mean(coal.size.unscaled, na.rm = T)), freq = F, col = "red", add = T)

# Gamma
lines(dgamma(min(coal.size.unscaled, na.rm = T):max(coal.size.unscaled, na.rm = T), shape = 3, scale = 1))
hist(round(rgamma(10000, 2, scale = 1)) + 1, freq = F, col = "blue", add = T)

# Negative binomial
lines(dnbinom(min(coal.size.unscaled, na.rm = T):max(coal.size.unscaled, na.rm = T), 30, mu = 1))
hist(rnbinom(10000, 8, mu = mean(coal.size.unscaled, na.rm = T)) + 1, freq = F, col = "green", add = T)

# Normal
lines(dnorm(min(coal.size.unscaled, na.rm = T):max(coal.size.unscaled, na.rm = T), mean = 2, sd = 0.7))

coal.size.lambda = exp(poisson.model$coefficients)


# Number of nomadic coalitions in the home range of a pride or 
# a resident male coalition
nb.nm.coal.hr.unscaled = read.csv("Data/027_Covariate_NbNMCoalHR.csv", row.names = 1)
nb.nm.coal.hr.unscaled = as.matrix(nb.nm.coal.hr.unscaled)
range(nb.nm.coal.hr.unscaled, na.rm = T)
nb.nm.coal.hr = (nb.nm.coal.hr.unscaled - mean(nb.nm.coal.hr.unscaled, na.rm = T)) / 
  (2 * sd(nb.nm.coal.hr.unscaled, na.rm = T)) # Standardize covariate

# Distribution
summary(glm(c(nb.nm.coal.hr.unscaled) ~ 1, "poisson"))
summary(MASS::glm.nb(c(nb.nm.coal.hr.unscaled) ~ 1))
negbin.model = MASS::glm.nb(c(nb.nm.coal.hr.unscaled) ~ 1)

hist(nb.nm.coal.hr.unscaled, freq = F)

# Poisson
lines(dpois(min(nb.nm.coal.hr.unscaled, na.rm = T):max(nb.nm.coal.hr.unscaled, na.rm = T), 1))
hist(rpois(10000, 2.3), freq = F, col = "red", add = T)

# Gamma
lines(dgamma(min(nb.nm.coal.hr.unscaled, na.rm = T):max(nb.nm.coal.hr.unscaled, na.rm = T), shape = 1, scale = 2.5))
hist(round(rgamma(10000, shape = 1, scale = 2.2)), freq = F, col = "blue", add = T, breaks = 15)

# Negative binomial
lines(dnbinom(min(nb.nm.coal.hr.unscaled, na.rm = T):max(nb.nm.coal.hr.unscaled, na.rm = T), 3, mu = 1))
hist(rnbinom(100000, negbin_model$theta, mu = mean(nb.nm.coal.hr.unscaled, na.rm = T)), freq = F, col = "green", add = T)
hist(rnbinom(100000, negbin_model$theta, prob = 0.35), freq = F, col = "green", add = T)

nb.nm.coal.hr.theta = negbin.model$theta


## 1.5. Lions groups across time ----
# ------------------------------

lions.groups = read.csv("Data/03_LionsGroups.csv", row.names = 1)



###########################################################################
#
# 2. Formatting capture histories ----
#
###########################################################################

## 2.1. Getting the first and last capture occasion of each lion ----
# --------------------------------------------------------------

get_first = function(x) min(which(x != 13)) # First observed state
get_last  = function(x){
  
  if(any(x == 11)){which(x == 11)} # Dead recovery occasion
  else{length(x)} # Otherwise, last capture occasion
  
}

lions_first = apply(lions_ch, 1, get_first)
lions_last = apply(lions_ch, 1, get_last)


## 2.2. Subset data if needed ----
# ---------------------------

reduced_data = F
subset_size = 100

if(reduced_data){
  
  # Subset capture histories
  lions_ch = lions_ch[1:subset_size, ]
  
  # Subset individual covariates
  age = age[1:subset_size, ]
  age_unscaled = age_unscaled[1:subset_size, ]
  habitat = habitat[1:subset_size, ]

  # Get new first and last sighting
  lions_first = apply(lions_ch, 1, get_first)
  lions_last  = apply(lions_ch, 1, get_last)
}




###########################################################################
#
# 3. Preparing and fitting the model ----
#
###########################################################################

## 3.1. Preparing the model ----
# -------------------------

# Model code

lions_code_lionKFDHMM = nimbleCode({
  # -------------------------------------------------
  #
  # Parameters:
  # mu.s.sa1:     survival probability of a Subadult 1 
  # mu.s.sa2f:    survival probability of a Subadult 2 Female
  # mu.s.sa2m:    survival probability of a Subadult 2 Male
  # mu.s.af:      survival probability of an Adult Female 
  # mu.s.ym:      survival probability of a Young Male 
  # mu.s.nm:      survival probability of a Nomadic Male 
  # mu.s.rm:      survival probability of a Resident Male 
  
  # mu.emig.ym:   emigration probability of a Young Male
  # mu.t.ym.nm:   transition probability YM-NM
  # mu.takeover:  takeover probability of a Nomadic Male
  # mu.eviction:  eviction probability of a Resident Male

  # mu.dp.pride:  detection probability of a pride individual (SA, AF, YM, RM)
  # mu.dp.nm:     detection probability of a Nomadic Male
  #
  # -------------------------------------------------
  #
  # States (S):
  #
  # 1 Subadult 1
  # 2 Subadult 2 Female
  # 3 Subadult 2 Male
  # 4 Adult Female
  # 5 Young Male 1
  # 6 Young Male 2
  # 7 Young Male 3
  # 8 Young Male 4
  # 9 Nomadic Male
  # 10 Resident Male
  # 11 Newly dead
  # 12 Permanently dead
  #
  # Observations (O): 
  #
  # 1 seen as Subadult 1
  # 2 seen as Subadult 2 Female
  # 3 seen as Subadult 2 Male
  # 4 seen as Adult Female
  # 5 seen as Young Male 1
  # 6 seen as Young Male 2
  # 7 seen as Young Male 3
  # 8 seen as Young Male 4
  # 9 seen as Nomadic Male
  # 10 seen as Resident Male
  # 11 seen Dead
  # 13 not seen
  #
  # -------------------------------------------------
  #
  # Generalised linear models
  
  for(i in 1:nind){
    
    for(t in lions.first[i]:lions.last[i]){ # From the first to the last sighting
                                            # of each individual
      
      # Survival
      
      # Young subadult
      logit(survSA1[i, t]) <- mu.s.sa1[season[t]] +
                              s.sa1.beta.nb.nm.coal.hr[season[t]] * nb.nm.coal.hr[group[i, t], t] +
                              s.sa1.beta.nb.af.pride[season[t]] * nb.af.pride[group[i, t], t] +
                              s.sa1.beta.habitat.woodland[season[t], habitat[i, t]] +
                              epsilon.s.sa1[season[t], year[t]]
      
      # Female old subadult
      logit(survSA2F[i, t]) <- mu.s.sa2f[season[t]] +
                               s.sa2.beta.nb.nm.coal.hr * nb.nm.coal.hr[group[i, t], t] +
                               s.sa2.beta.nb.af.pride[season[t]] * nb.af.pride[group[i, t], t] +
                               s.sa2.beta.habitat.woodland[season[t], habitat[i, t]] +
                               epsilon.s.sa2[season[t], year[t]]
      
      # Male old subadult                 
      logit(survSA2M[i, t]) <- mu.s.sa2m[season[t]] +
                               s.sa2.beta.nb.nm.coal.hr * nb.nm.coal.hr[group[i, t], t] +
                               s.sa2.beta.nb.af.pride[season[t]] * nb.af.pride[group[i, t], t] +
                               s.sa2.beta.habitat.woodland[season[t], habitat[i, t]] +
                               epsilon.s.sa2[season[t], year[t]]
      
      # Adult female
      logit(survAF[i, t]) <- mu.s.af[season[t]] +
                             s.af.beta.age[season[t]] * age[i, t] +
                             s.af.beta.nb.af.pride[season[t]] * nb.af.pride[group[i, t], t] +
                             s.af.beta.nb.nm.coal.hr[season[t]] * nb.nm.coal.hr[group[i, t], t] +
                             s.af.beta.habitat.woodland[season[t], habitat[i, t]] +
                             epsilon.s.af[season[t], year[t]]
      
      # Young male
      logit(survYM[i, t]) <- mu.s.ym[season[t]] +
                             s.ym.beta.nb.nm.coal.hr[season[t]] * nb.nm.coal.hr[group[i, t], t] +
                             s.ym.beta.nb.af.pride[season[t]] * nb.af.pride[group[i, t], t] +
                             s.ym.beta.habitat.woodland[season[t], habitat[i, t]] +
                             epsilon.s.ym[season[t], year[t]]
      
      # Nomadic male
      logit(survNM[i, t]) <- mu.s.nm[season[t]] +
                             s.nm.beta.coal.size[season[t]] * coal.size[group[i, t], t] +
                             s.nm.beta.habitat.woodland[season[t], habitat[i, t]] +
                             epsilon.s.nm[season[t], year[t]]
      
      # Resident male
      logit(survRM[i, t]) <- mu.s.rm[season[t]] +
                             s.rm.beta.coal.size[season[t]] * coal.size[group[i, t], t] +
                             s.rm.beta.nb.nm.coal.hr[season[t]] * nb.nm.coal.hr[group[i, t], t] +
                             s.rm.beta.habitat.woodland[season[t], habitat[i, t]] +
                             epsilon.s.rm[season[t], year[t]]
      
      
      # Transitions
      
      # Young-male emigration
      logit(emigYM[i, t]) <- mu.emig.ym[season[t]] +
                             emig.ym.beta.habitat.woodland[season[t], habitat[i, t]] +
                             epsilon.emig.ym[season[t], year[t]]
      
      # Young-male transition to nomadic male
      logit(transYMNM[i, t]) <- mu.t.ym.nm[season[t]] +
                                epsilon.t.ym.nm[season[t], year[t]]
      
      # Nomadic-male takeover
      logit(takeover[i, t]) <- mu.takeover[season[t]] +
                               takeover.beta.coal.size[season[t]] * coal.size[group[i, t], t] +
                               takeover.beta.habitat.woodland[season[t], habitat[i, t]] +
                               epsilon.takeover[season[t], year[t]]
      
      # Resident-male eviction
      logit(eviction[i, t]) <- mu.eviction[season[t]] +
                               eviction.beta.coal.size[season[t]] * coal.size[group[i, t], t] +
                               eviction.beta.nb.nm.coal.hr[season[t]] * nb.nm.coal.hr[group[i, t], t] +
                               eviction.beta.habitat.woodland[season[t], habitat[i, t]] +
                               epsilon.eviction[season[t], year[t]]
      
      
      # Detection probabilities
      
      # Pride individuals (subadults, young males, adult females, and resident males)
      logit(dpPride[i, t]) <- mu.dp.pride[season[t]] +
                              dp.pride.beta.habitat.woodland[season[t], habitat[i, t]] +
                              epsilon.dp.pride[season[t], year[t]]
      
      # Nomadic males
      logit(dpNM[i, t]) <- mu.dp.nm[season[t]] +
                           dp.nm.beta.habitat.woodland[season[t], habitat[i, t]] +
                           epsilon.dp.nm[season[t], year[t]]
      
      # Dead recovery
      logit(dpDead[i, t]) <- mu.dp.dead[season[t]] +
                             epsilon.dp.dead[season[t], year[t]]
      
    }
  }
  
  
  # Priors and constraints
  
  s.sa2.beta.nb.nm.coal.hr ~ dunif(-10, 10) # We do not estimate season-specific
                                            # effects of nomadic coalitions on
                                            # old subadults because of a lack
                                            # of data.

  # For the intercepts and all the other effects, we estimate one parameter 
  # per season
  for(seas in 1:2){
    
    # Means
    
    # Survival
    
    mu.s.sa1[seas] ~ dunif(-10, 10)     # Young subadult
    mu.s.sa2f[seas] ~ dunif(-10, 10)    # Female old subadult
    mu.s.sa2m[seas] ~ dunif(-10, 10)    # Male old subadult
    mu.s.af[seas] ~ dunif(-10, 10)      # Adult female
    mu.s.ym[seas] ~ dunif(-10, 10)      # Young male
    mu.s.nm[seas] ~ dunif(-10, 10)      # Nomadic male
    mu.s.rm[seas] ~ dunif(-10, 10)      # Resident male
    
    
    # Transition
    
    mu.emig.ym[seas] ~ dunif(-10, 10)   # Young-male emigration
    mu.t.ym.nm[seas] ~ dunif(-10, 10)   # Young-male transition to nomadic male
    mu.takeover[seas] ~ dunif(-10, 10)  # Nomadic-male takeover
    mu.eviction[seas] ~ dunif(-10, 10)  # Resident-male eviction
    
    
    # Detection probability
    
    mu.dp.pride[seas] ~ dunif(-10, 10)  # Pride individual
    mu.dp.nm[seas] ~ dunif(-10, 10)     # Nomadic male
    mu.dp.dead[seas] ~ dunif(-10, 10)   # Dead recovery
    
    
    # Other parameters
    
    # Young-subadult (SA1) survival
    
    s.sa1.beta.nb.nm.coal.hr[seas] ~ dunif(-10, 10)          # Number of nomadic coalitions in the home range
    s.sa1.beta.nb.af.pride[seas] ~ dunif(-10, 10)            # Number of adult females in the pride
    s.sa1.beta.habitat.woodland[seas, 1] <- 0                # Habitat - grassland (reference level)
    s.sa1.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)    # Habitat - woodland 

    
    # Old-subadult (SA2F and SA2M) survival
    # For the covariate effects we do not differentiate between sexes,
    # we do so only for the intercepts.
    
    s.sa2.beta.nb.af.pride[seas] ~ dunif(-10, 10)            # Number of adult females in the pride
    s.sa2.beta.habitat.woodland[seas, 1] <- 0                # Habitat - grassland (reference level)
    s.sa2.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)    # Habitat - woodland 

    
    # Adult-female (AF) survival
    
    s.af.beta.age[seas] ~ dunif(-10, 10)                     # Age
    s.af.beta.nb.af.pride[seas] ~ dunif(-10, 10)             # Number of adult females in the pride
    s.af.beta.nb.nm.coal.hr[seas] ~ dunif(-10, 10)           # Number of nomadic coalitions in the home range
    s.af.beta.habitat.woodland[seas, 1] <- 0                 # Habitat - grassland (reference level)
    s.af.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)     # Habitat - woodland 

    
    # Young-male (YM) survival
    
    s.ym.beta.nb.nm.coal.hr[seas] ~ dunif(-10, 10)           # Number of nomadic coalitions in the home range
    s.ym.beta.nb.af.pride[seas] ~ dunif(-10, 10)             # Number of adult females in the pride
    s.ym.beta.habitat.woodland[seas, 1] <- 0                 # Habitat - grassland (reference level)
    s.ym.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)     # Habitat - woodland 

    
    # Nomadic-male (NM) survival
    
    s.nm.beta.coal.size[seas] ~ dunif(-10, 10)               # Coalition size
    s.nm.beta.habitat.woodland[seas, 1] <- 0                 # Habitat - grassland (reference level)
    s.nm.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)     # Habitat - woodland 

    
    # Resident-male (RM) survival
    
    s.rm.beta.coal.size[seas] ~ dunif(-10, 10)               # Coalition size
    s.rm.beta.nb.nm.coal.hr[seas] ~ dunif(-10, 10)           # Number of nomadic coalitions in the home range
    s.rm.beta.habitat.woodland[seas, 1] <- 0                 # Habitat - grassland (reference level)
    s.rm.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)     # Habitat - woodland

    
    # Young-male (YM) emigration
    
    emig.ym.beta.habitat.woodland[seas, 1] <- 0              # Habitat - grassland (reference level)
    emig.ym.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)  # Habitat - woodland

    
    # Nomadic-male (NM) takeover
    
    takeover.beta.coal.size[seas] ~ dunif(-10, 10)           # Coalition size
    takeover.beta.habitat.woodland[seas, 1] <- 0             # Habitat - grassland (reference level)
    takeover.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10) # Habitat - woodland

    
    # Resident-male (RM) eviction
    
    eviction.beta.coal.size[seas] ~ dunif(-10, 10)           # Coalition size
    eviction.beta.nb.nm.coal.hr[seas] ~ dunif(-10, 10)       # Number of nomadic coalitions in the home range
    eviction.beta.habitat.woodland[seas, 1] <- 0             # Habitat - grassland (reference level) 
    eviction.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10) # Habitat - woodland
    
    
    # Detection probabilities
    
    # Pride individuals
    dp.pride.beta.habitat.woodland[seas, 1] <- 0             # Habitat - grassland (reference level) 
    dp.pride.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10) # Habitat - woodland
    
    # Nomadic males
    dp.nm.beta.habitat.woodland[seas, 1] <- 0                # Habitat - grassland (reference level) 
    dp.nm.beta.habitat.woodland[seas, 2] ~ dunif(-10, 10)    # Habitat - woodland
    
    # We only estimated season-specific intercepts for the detection probability
    # of dead individuals.
  }
  
  
  # Season-specific yearly random effects
  
  for(seas in 1:2){
    
    # Survival
    sigma.s.sa1[seas] ~ dunif(0, 10)     # Young subadult
    sigma.s.sa2[seas] ~ dunif(0, 10)     # Old subadult (both female and male)
    sigma.s.af[seas] ~ dunif(0, 10)      # Adult female
    sigma.s.ym[seas] ~ dunif(0, 10)      # Young male
    sigma.s.nm[seas] ~ dunif(0, 10)      # Nomadic male
    sigma.s.rm[seas] ~ dunif(0, 10)      # Resident male
    
    # Transitions
    sigma.emig.ym[seas] ~ dunif(0, 10)   # Young-male emigration
    sigma.t.ym.nm[seas] ~ dunif(0, 10)   # Young-male transition to nomadic male
    sigma.takeover[seas] ~ dunif(0, 10)  # Nomadic-male takeover
    sigma.eviction[seas] ~ dunif(0, 10)  # Resident-male eviction
    
    # Detection probabilities
    sigma.dp.pride[seas] ~ dunif(0, 10)  # Pride individual
    sigma.dp.nm[seas] ~ dunif(0, 10)     # Nomadic male
    sigma.dp.dead[seas] ~ dunif(0, 10)   # Dead recovery
    
  }
  
  for(yr in 1:30){
    
    for(seas in 1:2){
      
      # Survival
      epsilon.s.sa1[seas, yr] ~ dnorm(0, sd = sigma.s.sa1[seas])       # Young subadult
      epsilon.s.sa2[seas, yr] ~ dnorm(0, sd = sigma.s.sa2[seas])       # Old subadult (both female and male)
      epsilon.s.af[seas, yr] ~ dnorm(0, sd = sigma.s.af[seas])         # Adult female
      epsilon.s.ym[seas, yr] ~ dnorm(0, sd = sigma.s.ym[seas])         # Young male
      epsilon.s.nm[seas, yr] ~ dnorm(0, sd = sigma.s.nm[seas])         # Nomadic male
      epsilon.s.rm[seas, yr] ~ dnorm(0, sd = sigma.s.rm[seas])         # Resident male
      
      # Transitions
      epsilon.emig.ym[seas, yr] ~ dnorm(0, sd = sigma.emig.ym[seas])   # Young-male emigration
      epsilon.t.ym.nm[seas, yr] ~ dnorm(0, sd = sigma.t.ym.nm[seas])   # Young-male transition to nomadic male
      epsilon.takeover[seas, yr] ~ dnorm(0, sd = sigma.takeover[seas]) # Nomadic-male takeover
      epsilon.eviction[seas, yr] ~ dnorm(0, sd = sigma.eviction[seas]) # Resident-male eviction
      
      # Detection probabilities
      epsilon.dp.pride[seas, yr] ~ dnorm(0, sd = sigma.dp.pride[seas]) # Pride individual
      epsilon.dp.nm[seas, yr] ~ dnorm(0, sd = sigma.dp.nm[seas])       # Nomadic male
      epsilon.dp.dead[seas, yr] ~ dnorm(0, sd = sigma.dp.dead[seas])   # Dead recovery
      
    }
  }
  
  
  # Sampling of missing covariates values using the observed distributions
  
  for(t in 1:60){
    
    # Group-specific covariates
    for(g in 1:n_groups){

      # Number of females in the pride
      unscaled.nb.af.pride[g, t] ~ T(dnegbin(prob = 0.65, 
                                             size = nb.af.pride.theta), 
                                     min.nb.af.pride, max.nb.af.pride)
      nb.af.pride[g, t] <- (unscaled.nb.af.pride[g, t] - mu.nb.af.pride) / 
                           (2 * sd.nb.af.pride) # Standardize covariate

      # Coalition size
      unscaled.coal.size[g, t] ~ T(dpois(coal.size.lambda), 
                                   min.coal.size, max.coal.size)
      coal.size[g, t] <- (unscaled.coal.size[g, t] - mu.coal.size) / 
                         (2 * sd.coal.size) # Standardize covariate

      # Number of nomadic coalitions in the home range
      unscaled.nb.nm.coal.hr[g, t] ~ T(dnegbin(prob = 0.35, 
                                               size = nb.nm.coal.hr.theta), 
                                       min.nb.nm.coal.hr, max.nb.nm.coal.hr)
      nb.nm.coal.hr[g, t] <- (unscaled.nb.nm.coal.hr[g, t] - mu.nb.nm.coal.hr) / 
                             (2 * sd.nb.nm.coal.hr) # Standardize covariate

    }
    
    # Individual-specific covariate: Habitat
    for(i in 1:nind){
      
      habitat[i, t] <- temp.habitat[i, t] + 1
      temp.habitat[i, t] ~ dbin(size = 1, prob = habitat.prob)
      
    }
  }
  
  
  # Likelihood 
  
  for (i in 1:nind){
    for(u in 1:11){
      init[i, u] <- y2[i, lions.first[i]] == u # Get initial state distribution
    }
  }
  
  # Sample from custom distribution to assign state from first to last lion sighting
  for (i in 1:nind){
    y[i, lions.first[i]:lions.last[i]] ~ dDHMM_lionKF(length = (lions.last[i] - lions.first[i] + 1),            # Length of states vector
                                                      init = init[i, 1:11],                                     # Initial state distribution
                                                      survSA1 = survSA1[i, lions.first[i]:(lions.last[i])],     # Young-subadult survival
                                                      survSA2F = survSA2F[i, lions.first[i]:(lions.last[i])],   # Female old-subadult survival
                                                      survSA2M = survSA2M[i, lions.first[i]:(lions.last[i])],   # Male old-subadult survival
                                                      survAF = survAF[i, lions.first[i]:(lions.last[i])],       # Adult-female survival
                                                      survYM = survYM[i, lions.first[i]:(lions.last[i])],       # Young-male survival
                                                      survNM = survNM[i, lions.first[i]:(lions.last[i])],       # Nomadic-male survival
                                                      survRM = survRM[i, lions.first[i]:(lions.last[i])],       # Resident-male survival
                                                      emigYM = emigYM[i, lions.first[i]:(lions.last[i])],       # Nomadic-male emigration
                                                      transYMNM = transYMNM[i, lions.first[i]:(lions.last[i])], # Young-male transition to nomadic male
                                                      takeover = takeover[i, lions.first[i]:(lions.last[i])],   # Nomadic-male takeover
                                                      eviction = eviction[i, lions.first[i]:(lions.last[i])],   # Resident-male eviction
                                                      dpPride = dpPride[i, lions.first[i]:lions.last[i]],       # Detection of pride individuals
                                                      dpNM = dpNM[i, lions.first[i]:lions.last[i]],             # Detection of nomadic males
                                                      dpDead = dpDead[i, lions.first[i]:lions.last[i]])}        # Detection of dead individuals
})


# Model constants

lions_constants_lionKFDHMM = list(lions.first = lions.first,                                  # First sighting of each individual
                                   lions.last = lions.last,                                    # Last sighting of each individual
                                   n.occasions = ncol(lions.ch),                               # Total number of capture occasions
                                   nind = nrow(lions.ch),                                      # Total number of individuals
                                   n_groups = nrow(nb.af.pride.unscaled),                      # Number of lion groups
                                   group = apply(lions.groups,                                 # Group to which an individual belongs in a given timestep
                                                 c(1, 2), 
                                                 FUN = function(x) {
                                                   ifelse(!is.na(x), 
                                                          which(rownames(nb.af.pride.unscaled) == x), 
                                                          NA) 
                                                 }),   
                                   year = year,                                                # Vector of years
                                   season = season,                                            # Vector of seasons
                                   age = age,                                                  # Age of each individual at each timestep
                                   
                                   # Parameters for observed distributions and standarization
                                   nb.af.pride.theta = nb.af.pride.theta,                      # Theta for observed negative binomial distribution
                                                                                               # of the number of females in a pride
                                   min.nb.af.pride = min(nb.af.pride.unscaled, na.rm = T),     # Minimum number of females in a pride
                                                                                               # for observed negative binomial distribution
                                   max.nb.af.pride = max(nb.af.pride.unscaled, na.rm = T),     # Maximum number of females in a pride
                                                                                               # for observed negative binomial distribution
                                   mu.nb.af.pride = mean(nb.af.pride.unscaled, na.rm = T),     # Mean number of females in a pride for standardization
                                   sd.nb.af.pride = sd(nb.af.pride.unscaled, na.rm = T),       # Standard deviation of number of females in a pride for standardization
                                   coal.size.lambda = coal.size.lambda,                        # Lambda for observed Poisson distribution of coalition size
                                   min.coal.size = min(coal.size.unscaled, na.rm = T),         # Minimun coalition size for observed Poisson distribution
                                   max.coal.size = max(coal.size.unscaled, na.rm = T),         # Maximun coalition size for observed Poisson distribution
                                   mu.coal.size = mean(coal.size.unscaled, na.rm = T),         # Mean coalition size for standardization
                                   sd.coal.size = sd(coal.size.unscaled, na.rm = T),           # Standard deviation of coalition size for standardization
                                   nb.nm.coal.hr.theta = nb.nm.coal.hr.theta,                  # Theta for observed negative binomial distribution
                                                                                               # of the number of nomadic coalitions in the home range 
                                   min.nb.nm.coal.hr = min(nb.nm.coal.hr.unscaled, na.rm = T), # Minimum number of nomadic coalitions in the home range
                                                                                               # for observed negative binomial distribution
                                   max.nb.nm.coal.hr = max(nb.nm.coal.hr.unscaled, na.rm = T), # Maximum number of nomadic coalitions in the home range
                                                                                               # for observed negative binomial distribution
                                   mu.nb.nm.coal.hr = mean(nb.nm.coal.hr.unscaled, na.rm = T), # Mean number of nomadic coalitions in the home range
                                                                                               # for standardization
                                   sd.nb.nm.coal.hr = sd(nb.nm.coal.hr.unscaled, na.rm = T),   # Standard deviation of the number of nomadic 
                                                                                               # coalitions in the home range for standardization
                                   habitat.prob = inv.logit(habitat.intercept.estimate))       # Observed probability of being in the woodland habitat


# Model data

lions_data_lionKFDHMM = list(y = lions.ch,                                    # Capture histories
                              y2 = lions.ch,                                   # Capture histories for initial state distribution
                              unscaled.nb.af.pride = nb.af.pride.unscaled,     # Non-standardized number of adult females in a pride 
                                                                               # for the sampling of missing covariate values
                              nb.af.pride = nb.af.pride,                       # Standardized number of adult females in a pride for the model
                              unscaled.coal.size = coal.size.unscaled,         # Non-standardized coalition size for the sampling 
                                                                               # of missing covariate values
                              coal.size = coal.size,                           # Standardized coalition size for the model
                              unscaled.nb.nm.coal.hr = nb.nm.coal.hr.unscaled, # Non-standardized number of nomadic coalitions in 
                                                                               # the home range for the sampling of missing covariate values
                              nb.nm.coal.hr = nb.nm.coal.hr,                   # Standardized number of nomadic coalitions in the home range for the model
                              temp.habitat = habitat - 1,                      # Habitat as a (0, 1) variable for the sampling of missing covariate values
                              habitat = habitat)                               # Habitat as a (1, 2) variable for the model


# Model initial values

# Initial values for the sampling of covariates with missing values

# Number of adult females in a pride
inits.nb.af.pride = function(cov){ # cov is the observed data
  
  cov.new = matrix(NA, nrow = nrow(cov), ncol = ncol(cov)) # Empty matrix
  
  # Replace NAs by sampling the observed negative binomial distribution
  for(i in 1:nrow(cov.new)){
    
    which.na = which(is.na(cov[i, ]))
    if(length(which.na) == 0){next}
    cov.new[i, which.na] = rnbinom(length(which.na), size = nb.af.pride.theta, prob = 0.65)
    
  }
  
  # Bound the sampled values to the minimum and maximum observed values
  cov.new[which(cov.new > max(cov, na.rm = T))] = max(cov, na.rm = T)
  cov.new[which(cov.new < min(cov, na.rm = T))] = min(cov, na.rm = T)
  
  return(cov.new) # Full covariate data
}


# Coalition size
inits.coal.size = function(cov){ # cov is the observed data
  
  cov.new = matrix(NA, nrow = nrow(cov), ncol = ncol(cov)) # Empty matrix
  
  # Replace NAs by sampling the observed Poisson distribution
  for(i in 1:nrow(cov.new)){
    
    which.na = which(is.na(cov[i, ]))
    if(length(which.na) == 0){next}
    cov.new[i, which.na] = rpois(length(which.na), coal.size.lambda) + 1
    
  }
  
  # Bound the sampled values to the minimum and maximum observed values
  cov.new[which(cov.new > max(cov, na.rm = T))] = max(cov, na.rm = T)
  cov.new[which(cov.new < min(cov, na.rm = T))] = min(cov, na.rm = T)
  
  return(cov.new) # Full covariate data
}


# Number of nomadic coalitions in the home range
inits.nb.nm.coal.hr = function(cov){ # cov is the observed data
  
  cov.new = matrix(NA, nrow = nrow(cov), ncol = ncol(cov)) # Empty matrix
  
  # Replace NAs by sampling the observed negative binomial distribution
  for(i in 1:nrow(cov.new)){
    
    which.na = which(is.na(cov[i, ]))
    if(length(which.na) == 0){next}
    cov.new[i, which.na] = rnbinom(length(which.na), size = nb.nm.coal.hr.theta, prob = 0.35)
    
  }
  
  # Bound the sampled values to the minimum and maximum observed values
  cov.new[which(cov.new > max(cov, na.rm = T))] = max(cov, na.rm = T)
  cov.new[which(cov.new < min(cov, na.rm = T))] = min(cov, na.rm = T)
  
  return(cov.new) # Full covariate data
}


# Habitat 
inits.habitat = function(cov){ # cov is the observed data
  
  cov.new = matrix(NA, nrow = nrow(cov), ncol = ncol(cov)) # Empty matrix
  
  # Replace NAs by sampling the observed binomial distribution
  for(i in 1:nrow(cov.new)){
    
    which.na = which(is.na(cov[i, ]))
    if(length(which.na) == 0){next}
    cov.new[i, which.na] = rbinom(length(which.na), size = 1, prob = 0.3)
    
  }  
  
  return(cov.new) # Full covariate data
}


# Setting initial values for the model parameters and creating the list
# of initial values

lions_inits_lionKFDHMM = function(){
  list(# Intercepts
       mu.s.sa1 = runif(2, -2, 2),      # Young-subadult survival
       mu.s.sa2f = runif(2, -2, 2),     # Female old-subadult survival 
       mu.s.sa2m = runif(2, -2, 2),     # Male old-subadult survival
       mu.s.af = runif(2, -2, 2),       # Adult-female survival
       mu.s.ym = runif(2, -2, 2),       # Young-male survival
       mu.s.nm = runif(2, -2, 2),       # Nomadic-male survival
       mu.s.rm = runif(2, -2, 2),       # Resident-male survival
       mu.emig.ym = runif(2, -2, 2),    # Young-male emigration
       mu.t.ym.nm = runif(2, -2, 2),    # Young-male transition to nomadic male
       mu.takeover = runif(2, -2, 2),   # Nomadic-male takeover
       mu.eviction = runif(2, -2, 2),   # Resident-male eviction 
       mu.dp.pride = runif(2, -2, 2),   # Detection of pride individual
       mu.dp.nm = runif(2, -2, 2),      # Detection of nomadic male
       mu.dp.dead = runif(2, -2, 2),    # Detection of dead individual
       
       # Betas
       # Young-subadult survival 
       s.sa1.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),       # Number of nomadic coalitions in the home range
       s.sa1.beta.nb.af.pride = runif(2, -0.05, 0.05),         # Number of adult females in the pride
       s.sa1.beta.habitat.woodland = matrix(c(rep(NA, 2),      # Habitat (NA for grassland which is the reference level)
                                              runif(2, 
                                                    -0.05, 
                                                    0.05)), 
                                            2, 2),
       # Old-subadult survival 
       s.sa2.beta.nb.nm.coal.hr = runif(1, -0.05, 0.05),       # Number of nomadic coalitions in the home range 
       s.sa2.beta.nb.af.pride = runif(2, -0.05, 0.05),         # Number of adult females in the pride
       s.sa2.beta.habitat.woodland = matrix(c(rep(NA, 2),      # Habitat (NA for grassland which is the reference level)
                                              runif(2, 
                                                    -0.05, 
                                                    0.05)), 
                                            2, 2),
       # Adult-female survival
       s.af.beta.age = runif(2, -0.05, 0.05),                  # Age
       s.af.beta.nb.af.pride = runif(2, -0.05, 0.05),          # Number of adult females in the pride
       s.af.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),        # Number of nomadic coalitions in the home range
       s.af.beta.habitat.woodland = matrix(c(rep(NA, 2),       # Habitat (NA for grassland which is the reference level)
                                             runif(2, 
                                                   -0.05, 
                                                   0.05)), 
                                           2, 2),
       # Young-male survival
       s.ym.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),       # Number of nomadic coalitions in the home range
       s.ym.beta.nb.af.pride = runif(2, -0.05, 0.05),         # Number of adult females in the pride
       s.ym.beta.habitat.woodland = matrix(c(rep(NA, 2),      # Habitat (NA for grassland which is the reference level)
                                             runif(2, 
                                                   -0.05, 
                                                   0.05)), 
                                           2, 2),
       # Nomadic-male survival
       s.nm.beta.coal.size = runif(2, -0.05, 0.05),           # Coalition size
       s.nm.beta.habitat.woodland = matrix(c(rep(NA, 2),      # Habitat (NA for grassland which is the reference level)
                                             runif(2, 
                                                   -0.05, 
                                                   0.05)), 
                                           2, 2),
       s.rm.beta.coal.size = runif(2, -0.05, 0.05),           # Coalition size
       s.rm.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),       # Number of nomadic coalitions in the home range
       s.rm.beta.habitat.woodland = matrix(c(rep(NA, 2),      # Habitat (NA for grassland which is the reference level)
                                             runif(2, 
                                                   -0.05, 
                                                   0.05)), 
                                           2, 2),
       # Young-male emigration
       emig.ym.beta.habitat.woodland = matrix(c(rep(NA, 2),   # Habitat (NA for grassland which is the reference level) 
                                                runif(2, 
                                                      -0.05, 
                                                      0.05)), 
                                              2, 2),
       # Nomadic-male takeover
       takeover.beta.coal.size = runif(2, -0.05, 0.05),       # Coalition size
       takeover.beta.habitat.woodland = matrix(c(rep(NA, 2),  # Habitat (NA for grassland which is the reference level) 
                                                 runif(2, 
                                                       -0.05, 
                                                       0.05)), 
                                               2, 2),
       # Resident-male eviction
       eviction.beta.coal.size = runif(2, -0.05, 0.05),       # Coalition size
       eviction.beta.nb.nm.coal.hr = runif(2, -0.05, 0.05),   # Number of nomadic coalitions in the home range
       eviction.beta.habitat.woodland = matrix(c(rep(NA, 2),  # Habitat (NA for grassland which is the reference level) 
                                                 runif(2, 
                                                       -0.05, 
                                                       0.05)), 
                                               2, 2),
       # Pride-individual detection 
       dp.pride.beta.habitat.woodland = matrix(c(rep(NA, 2),  # Habitat (NA for grassland which is the reference level) 
                                                 runif(2, 
                                                       -0.05, 
                                                       0.05)), 
                                               2, 2),
       # Nomadic-male detection
       dp.nm.beta.habitat.woodland = matrix(c(rep(NA, 2),     # Habitat (NA for grassland which is the reference level) 
                                              runif(2, 
                                                    -0.05, 
                                                    0.05)), 
                                            2, 2),
       
       # Random effect sigmas
       sigma.s.sa1 = runif(2, 0.05, 0.5),    # Young-subadult survival
       sigma.s.sa2 = runif(2, 0.05, 0.5),    # Old-subadult survival
       sigma.s.af = runif(2, 0.05, 0.5),     # Adult-female survival
       sigma.s.ym = runif(2, 0.05, 0.5),     # Young-male survival
       sigma.s.nm = runif(2, 0.05, 0.5),     # Nomadic-male survival
       sigma.s.rm = runif(2, 0.05, 0.5),     # Resident-male survival
       sigma.emig.ym = runif(2, 0.05, 0.5),  # Young-male emigration
       sigma.t.ym.nm = runif(2, 0.05, 0.5),  # Young-male transition to nomadic male
       sigma.takeover = runif(2, 0.05, 0.5), # Nomadic-male takeover
       sigma.eviction = runif(2, 0.05, 0.5), # Resident-male eviction
       sigma.dp.pride = runif(2, 0.05, 0.5), # Pride-individual detection
       sigma.dp.nm = runif(2, 0.05, 0.5),    # Nomadic-male detection
       sigma.dp.dead = runif(2, 0.05, 0.5),  # Dead-individual detection
       
       # Random effect epsilons
       epsilon.s.sa1 = matrix(0, nrow = max(season), ncol = max(year)),    # Young-subadult survival  
       epsilon.s.sa2 = matrix(0, nrow = max(season), ncol = max(year)),    # Old-subadult survival
       epsilon.s.af = matrix(0, nrow = max(season), ncol = max(year)),     # Adult-female survival
       epsilon.s.ym = matrix(0, nrow = max(season), ncol = max(year)),     # Young-male survival
       epsilon.s.nm = matrix(0, nrow = max(season), ncol = max(year)),     # Nomadic-male survival
       epsilon.s.rm = matrix(0, nrow = max(season), ncol = max(year)),     # Resident-male survival
       epsilon.emig.ym = matrix(0, nrow = max(season), ncol = max(year)),  # Young-male emigration
       epsilon.t.ym.nm = matrix(0, nrow = max(season), ncol = max(year)),  # Young-male transition to nomadic male
       epsilon.takeover = matrix(0, nrow = max(season), ncol = max(year)), # Nomadic-male takeover
       epsilon.eviction = matrix(0, nrow = max(season), ncol = max(year)), # Resident-male eviction
       epsilon.dp.pride = matrix(0, nrow = max(season), ncol = max(year)), # Pride-individual detection
       epsilon.dp.nm = matrix(0, nrow = max(season), ncol = max(year)),    # Nomadic-male detection
       epsilon.dp.dead = matrix(0, nrow = max(season), ncol = max(year)),  # Dead-individual detection
       
       # Missing covariate initial values
       unscaled.nb.af.pride = inits.nb.af.pride(nb.af.pride.unscaled),       # Non-standardized number of adult females in a pride
       unscaled.coal.size = inits.coal.size(coal.size.unscaled),             # Non-standardized coalition size
       unscaled.nb.nm.coal.hr = inits.nb.nm.coal.hr(nb.nm.coal.hr.unscaled), # Non-standardized number of nomadic coalitions in the home range
       temp.habitat = inits.habitat(habitat))                                # Habitat with (0, 1) values
}


# Sample initial values
init.val = list(lions_inits_lionKFDHMM(), # Chain 1
                lions_inits_lionKFDHMM(), # Chain 2
                lions_inits_lionKFDHMM(), # Chain 3
                lions_inits_lionKFDHMM()) # Chain 4


# Parameters monitored
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



## 3.2. Fitting the model ----
# -----------------------

# Function to parallelize the chains
startnimbleMCMC = function(sim){ # sim = computer core
  
  if (sim > 1) {      # Give each core a random seed
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
  }
  
  # Run 5000 iterations of the model
  nimbleMCMC(code = lions_code_lionKFDHMM,
             constants = lions_constants_lionKFDHMM,
             data = lions_data_lionKFDHMM,
             monitors = parameters.lionKFDHMM,
             inits = lions_inits_lionKFDHMM,
             niter = 5000,
             nburnin = 0,
             nchains = 1,
             thin = 1,
             samplesAsCodaMCMC = TRUE)
  
}

ncpus = 4 # Number of cores


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run.

# Set up parallel environment
library(snowfall)

sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/MCMCProgress.txt") # Initialisation and progress output file
sfExport(list = c(ls(), ".Random.seed")) # Export current environment to each core
sfLibrary(nimble, warn.conflicts = FALSE) # Load nimble 
sfLibrary("snowfall", character.only = TRUE) # Load snowfall

# Results of each set of iterations are saved to a file
res.iter.set = sfClusterApplyLB(1:ncpus, startnimbleMCMC) # Run chains on one core each

# Setup filename based on data and time and save output
filename_temp = paste0(sub('\\..*', '', "Output/LionsFullMultistateModel_Output"), format(Sys.time(),'_%m%d_%H%M%S'))
filename = paste0(filename_temp,'.RData')

save(res.iter.set, file = filename)
sfStop() # Stop parallelization