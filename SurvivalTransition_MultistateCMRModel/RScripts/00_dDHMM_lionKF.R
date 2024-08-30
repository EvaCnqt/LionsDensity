############################################################################
#
# The aim of this script is to create a custom distribution to fit a 
# multistate model on lions capture histories.
#
# This custom distribution is based on Nater et al. (2020), Journal of Animal
# Ecology ( https://doi.org/10.1111/1365-2656.13269)
#
# Author: Eva Conquet
###########################################################################


###########################################################################
#
# 1. House keeping and loading libraries and data
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading libraries ----
# -----------------------

library(nimble)

nimbleOptions(disallow_multivariate_argument_expressions = F, enableWAIC = T)




###########################################################################
#
# 2. Creating the custom HMM distribution
#
###########################################################################

# States (S):

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
# 11 Dead
# 12 Permanently dead


# Observations (O): 

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
# 11 seen as Dead
# 13 not seen


dDHMM_lionKF = nimbleFunction(
  run = function(
    
    ## Argument type declarations
    
    x = double(1),           # Vector containing the observed capture history data
    length = double(),       # Length of the capture history
    init = double(1),        # Initial state probabilities
    survSA1 = double(1),     # State-specific survival
    survSA2F = double(1),
    survSA2M = double(1),
    survAF = double(1),
    survYM = double(1),
    survNM = double(1),
    survRM = double(1),
    transYMNM = double(1),   # Between-state transitions
    emigYM = double(1),
    takeover = double(1),
    eviction = double(1),
    dpPride = double(1),       # Detection probabilities
    dpNM = double(1),
    dpDead = double(1),
    log = double()){         # Logical argument specifying whether the log of the likelihood should be returned
    
    logL <- 0                # Initialise log-likelihood
    pi <- init               # Initialise state probabilities
    
    for(t in 1:length){      # Iterate over observations
      
      # x = "recorded as"
      # pi = probability of each latent state, conditioned on preceding observations
      # Zpi = probability of current observed capture, conditioned on each possible latent state
      
      Zpi <- pi  # Initialise Zpi with the values of pi to avoid assigning values        
                 # to Zpi when the observation probability of a given latent
                 # state in a given observed state is 1 (e.g. Zpi[12] when x[t] == 13,
                 # as permanently dead individuals will always be unobserved).
      
      # Detection probabilities 
      
      if(x[t] == 1){ 
        
        # We do not assign any value to Zpi[1] here because the latent state 1 "young subadults"
        # is the first state defined in our model. Therefore, in the capture histories,
        # observations are either (1) an NA if the first capture of an individual took place
        # when it was older than 1.5 years, or (2) a 1 if the first capture happened when
        # it was between 1 and 1.5 years old.
        
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 2){ 
        
        Zpi[1] <- 0
        Zpi[2] <- pi[2] * dpPride[t]
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 3){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- pi[3] * dpPride[t]
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 4){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- pi[4] * dpPride[t]
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 5){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- pi[5] * dpPride[t]
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 6){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- pi[6] * dpPride[t]
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 7){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- pi[7] * dpPride[t]
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 8){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- pi[8] * dpPride[t]
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 9){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- pi[9] * dpNM[t]
        Zpi[10] <- 0
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 10){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- pi[10] * dpPride[t]
        Zpi[11] <- 0
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 11){ 
        
        Zpi[1] <- 0
        Zpi[2] <- 0
        Zpi[3] <- 0
        Zpi[4] <- 0
        Zpi[5] <- 0
        Zpi[6] <- 0
        Zpi[7] <- 0
        Zpi[8] <- 0
        Zpi[9] <- 0
        Zpi[10] <- 0
        Zpi[11] <- pi[11] * dpDead[t]
        Zpi[12] <- 0
        
      }
      
      if(x[t] == 13){ 
        
        Zpi[1] <- 0
        Zpi[2] <- pi[2] * (1 - dpPride[t])
        Zpi[3] <- pi[3] * (1 - dpPride[t])
        Zpi[4] <- pi[4] * (1 - dpPride[t])
        Zpi[5] <- pi[5] * (1 - dpPride[t])
        Zpi[6] <- pi[6] * (1 - dpPride[t])
        Zpi[7] <- pi[7] * (1 - dpPride[t])
        Zpi[8] <- pi[8] * (1 - dpPride[t])
        Zpi[9] <- pi[9] * (1 - dpNM[t])
        Zpi[10] <- pi[10] * (1 - dpPride[t])
        Zpi[11] <- pi[11] * (1- dpDead[t])
        
        # We do not assign any value to Zpi[12] here because individuals in
        # the latent "permanently dead" state 12 are always unobserved (observed state 13).
        # The value of Zpi[12] is therefore the one it has been initialised with (pi[12])
        
      }
      
      sumZpi <- sum(Zpi)         # Log-likelihood contribution of given observed state x
      logL <- logL + log(sumZpi) # Overall log likelihood
      
      # Transition probabilities 
      
      if(t != length){
        
        pi[1] <- 0
        pi[2] <- Zpi[1] * survSA1[t] * 0.55
        pi[3] <- Zpi[1] * survSA1[t] * (1 - 0.55)
        pi[4] <- Zpi[2] * survSA2F[t] + 
                 Zpi[4] * survAF[t]
        pi[5] <- Zpi[3] * survSA2M[t]
        pi[6] <- Zpi[5] * survYM[t] * (1 - emigYM[t])
        pi[7] <- Zpi[6] * survYM[t] * (1 - emigYM[t])
        pi[8] <- Zpi[7] * survYM[t] * (1 - emigYM[t])
        pi[9] <- Zpi[5] * survYM[t] * emigYM[t] * transYMNM[t] + 
                 Zpi[6] * survYM[t] * emigYM[t] * transYMNM[t] + 
                 Zpi[7] * survYM[t] * emigYM[t] * transYMNM[t] + 
                 Zpi[8] * survYM[t] * transYMNM[t]  + 
                 Zpi[9] * survNM[t] * (1 - takeover[t]) + 
                 Zpi[10] * survRM[t] * eviction[t]
        pi[10] <- Zpi[5] * survYM[t] * emigYM[t] * (1 - transYMNM[t]) + 
                  Zpi[6] * survYM[t] * emigYM[t] * (1 - transYMNM[t]) + 
                  Zpi[7] * survYM[t] * emigYM[t] * (1 - transYMNM[t]) + 
                  Zpi[8] * survYM[t] * (1 - transYMNM[t]) + 
                  Zpi[9] * survNM[t] * takeover[t] + 
                  Zpi[10] * survRM[t] * (1 - eviction[t])
        pi[11] <- Zpi[1] * (1 - survSA1[t]) + 
                  Zpi[2] * (1 - survSA2F[t]) + 
                  Zpi[3] * (1 - survSA2M[t]) + 
                  Zpi[4] * (1 - survAF[t]) + 
                  Zpi[5] * (1 - survYM[t]) + 
                  Zpi[6] * (1 - survYM[t]) + 
                  Zpi[7] * (1 - survYM[t]) + 
                  Zpi[8] * (1 - survYM[t]) + 
                  Zpi[9] * (1 - survNM[t]) + 
                  Zpi[10] * (1 - survRM[t]) 
        pi[12] <- Zpi[11] + Zpi[12]
        
        pi <- pi / sumZpi  # Normalise
      }
    }
    
    returnType(double())
    
    if(log) return(logL) else return(exp(logL)) # Return log-likelihood
  }
)


rDHMM_lionKF = nimbleFunction(
  run = function(
    n = double(),
    length = double(),
    init = double(1),
    survSA1 = double(1),      # State-specific survival
    survSA2F = double(1),
    survSA2M = double(1),
    survAF = double(1),
    survYM = double(1),
    survNM = double(1),
    survRM = double(1),
    transYMNM = double(1),    # Between-state transitions
    emigYM = double(1),
    takeover = double(1),
    eviction = double(1),
    dpPride = double(1),        # Detection probabilities
    dpNM = double(1),
    dpDead = double(1)){
    
    x <- rep(1, length)
    
    returnType(double(1))
    return(x)
    
  }
)

registerDistributions(list(
  dDHMM_lionKF = list(
    BUGSdist = 'dDHMM_lionKF(length, init, survSA1, survSA2F, survSA2M, survAF, survYM, survNM, survRM, emigYM, transYMNM, takeover, eviction, dpPride, dpNM, dpDead)',
    types = c('value = double(1)', 'init = double(1)', 'survSA1 = double(1)', 'survSA2F = double(1)', 'survSA2M = double(1)', 'survAF = double(1)', 'survYM = double(1)', 'survNM = double(1)', 'survRM = double(1)', 'emigYM = double(1)', 'transYMNM = double(1)', 'takeover = double(1)', 'eviction = double(1)', 'dpPride = double(1)', 'dpNM = double(1)', 'dpDead = double(1)'),
    discrete = TRUE
  )
))
