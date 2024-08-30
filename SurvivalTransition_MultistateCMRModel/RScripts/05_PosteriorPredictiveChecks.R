############################################################################
#
# This script uses samples obtained from chains of an MCMC algorithm.
#
# The aim of this script is to perform posterior predictive checks using
# simulated capture histories from sampled parameter values.
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

# Individual capture histories
lions.ch = read.csv("Data/011_LionsCaptureHistories.csv", row.names = 1)
lions.ch = as.matrix(lions.ch)

# Simulated datasets
load("Output/Lions_MultistateModel_Simulated_Data.RData")
data_simulation_output = c(lions_multistate_simulated_data[[1]], # Merge simulated datasets
                           lions_multistate_simulated_data[[2]], # from each parallel core
                           lions_multistate_simulated_data[[3]],
                           lions_multistate_simulated_data[[4]], 
                           lions_multistate_simulated_data[[5]])




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

lions.first = apply(lions.ch, 1, get_first)
lions.last = apply(lions.ch, 1, get_last)




###########################################################################
#
# 3. Functions to calculate summary measures on datasets ----
#
###########################################################################

## 3.1. Total number of recaptures in any state or a specific state  ----
# -----------------------------------------------------------------

nb.recaptures = function(data,
                         state){
  
  # Initialize the counter of first captures in the focal state
  # to remove them from the total number of recaptures in that state.
  # This is because we initialized the simulated capture histories
  # with the state observed at the first capture.
  nb_first_captures_in_state = 0
  
  # Add 1 to the counter if a first capture correspond to the focal state
  for(i in 1:nrow(data)){
    
    if(data[i, lions.first[i]] %in% state) nb_first_captures_in_state = 
        nb_first_captures_in_state + 1
  }
  
  # Return the number of recaptures in the focal state removing the number
  # of first captures in that state
  return(length(which(data %in% state)) - nb_first_captures_in_state)
  
}


## 3.2. Number of recaptures in any state or a specific state at t + n  ----
# --------------------------------------------------------------------

nb.recaptures.tn = function(data,
                            state,
                            n){
  
  # Vector of number of recaptures at t+1 for each t
  nb_recaptures_tn = c()
  
  for(i in 1:(ncol(data) - n)){
    
    # Subset data to individuals first seen in current capture occasion (i)
    data_subset = data[which(lions.first == i), ]
    
    # Add number of recaptures in the focal state at t+n to the vector
    nb_recaptures_tn = c(nb_recaptures_tn, length(which(data_subset[, i + n] %in% state)))
  }
  
  return(nb_recaptures_tn)
  
}


## 3.3. Number of recaptures in any state or a specific state from a given state at t + n  ----
# ---------------------------------------------------------------------------------------

nb.recaptures.from.state.tn = function(data,
                                       state_start,
                                       state_end,
                                       n){
  
  # Vector of number of recaptures at t+1 for each t
  nb_recaptures_from_state_tn = c()
  
  for(i in 1:(ncol(data) - n)){
    
    # Subset data to individuals first seen in current capture occasion (i)
    data_subset = data[which(lions.first == i), ]
    
    # Add number of recaptures in the focal state at t+n to the vector
    nb_recaptures_from_state_tn = c(nb_recaptures_from_state_tn, 
                                    length(which(data_subset[, i] %in% state_start & 
                                                   data_subset[, i + n] %in% state_end)))
  }
  
  return(nb_recaptures_from_state_tn)
  
}




###########################################################################
#
# 4. Calculate summary measures on true and simulated datasets ----
#
###########################################################################

## 4.1. True data ----
# ---------------

## 4.1.1. Total number of recaptures ----
# ----------------------------------

total_nb_recaptures_true = nb.recaptures(lions.ch, seq(1, 10))


## 4.1.2. Total number of recaptures at t + 1 ----
# -------------------------------------------

total_nb_recaptures_t1_true = nb.recaptures.tn(lions.ch, seq(1, 10), 1)


## 4.1.3. Total number of recaptures at t + 2 ----
# -------------------------------------------

total_nb_recaptures_t2_true = nb.recaptures.tn(lions.ch, seq(1, 10), 2)


## 4.1.4. Total number of recaptures as female old subadults (2) ----
# --------------------------------------------------------------

total_nb_recaptures_sa2f_true = nb.recaptures(lions.ch, 2)


## 4.1.5. Number of recaptures as female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------

total_nb_recaptures_sa2f_t1_true = nb.recaptures.tn(lions.ch, 2, 1)


## 4.1.6. Total number of recaptures as male old subadults (3) ----
# ------------------------------------------------------------

total_nb_recaptures_sa2m_true = nb.recaptures(lions.ch, 3)


## 4.1.7. Number of recaptures as male old subadults (3) at t + 1 ----
# ---------------------------------------------------------------

total_nb_recaptures_sa2m_t1_true = nb.recaptures.tn(lions.ch, 3, 1)


## 4.1.8. Total number of recaptures as young male (5, 6, 7, 8) ----
# -------------------------------------------------------------

total_nb_recaptures_ym_true = nb.recaptures(lions.ch, seq(5, 8))


## 4.1.9. Number of recaptures as young male (5, 6, 7, 8) at t + 1 ----
# ----------------------------------------------------------------

total_nb_recaptures_ym_t1_true = nb.recaptures.tn(lions.ch, seq(5, 8), 1)


## 4.1.10. Number of recaptures as young male (5, 6, 7, 8) at t + 2 ----
# -----------------------------------------------------------------

total_nb_recaptures_ym_t2_true = nb.recaptures.tn(lions.ch, seq(5, 8), 2)


## 4.1.11. Total number of recaptures as young male 1 (5) ----
# -------------------------------------------------------

total_nb_recaptures_ym1_true = nb.recaptures(lions.ch, 5)


## 4.1.12. Number of recaptures as young male 1 (5) at t + 1 ----
# ----------------------------------------------------------

total_nb_recaptures_ym1_t1_true = nb.recaptures.tn(lions.ch, 5, 1)


## 4.1.13. Number of recaptures as young male 1 (5) at t + 2 ----
# ----------------------------------------------------------

total_nb_recaptures_ym1_t2_true = nb.recaptures.tn(lions.ch, 5, 2)


## 4.1.14. Number of recaptures as young male 1 (5) from male old subadult (3) at t + 1 ----
# -------------------------------------------------------------------------------------

total_nb_recaptures_ym1_from_sa2m_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                        3, 5, 1)


## 4.1.15. Number of recaptures as young male 1 (5) from young subadult (1) at t + 2 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_ym1_from_sa1_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       1, 5, 2)


## 4.1.16. Total number of recaptures as young male 2 (6) ----
# -------------------------------------------------------

total_nb_recaptures_ym2_true = nb.recaptures(lions.ch, 6)


## 4.1.17. Number of recaptures as young male 2 (6) at t + 1 ----
# ----------------------------------------------------------

total_nb_recaptures_ym2_t1_true = nb.recaptures.tn(lions.ch, 6, 1)


## 4.1.18. Number of recaptures as young male 2 (6) at t + 2 ----
# ----------------------------------------------------------

total_nb_recaptures_ym2_t2_true = nb.recaptures.tn(lions.ch, 6, 2)


## 4.1.19. Number of recaptures as young male 2 (6) from young male 1 (5) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_ym2_from_ym1_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       5, 6, 1)


## 4.1.20. Number of recaptures as young male 2 (6) from male old subadult (3) at t + 2 ----
# -------------------------------------------------------------------------------------

total_nb_recaptures_ym2_from_sa2m_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                        3, 6, 2)


## 4.1.21. Total number of recaptures as young male 3 (7) ----
# -------------------------------------------------------

total_nb_recaptures_ym3_true = nb.recaptures(lions.ch, 7)


## 4.1.22. Number of recaptures as young male 3 (7) at t + 1 ----
# ----------------------------------------------------------

total_nb_recaptures_ym3_t1_true = nb.recaptures.tn(lions.ch, 7, 1)


## 4.1.23. Number of recaptures as young male 3 (7) at t + 2 ----
# ----------------------------------------------------------

total_nb_recaptures_ym3_t2_true = nb.recaptures.tn(lions.ch, 7, 2)


## 4.1.24. Number of recaptures as young male 3 (7) from young male 2 (6) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_ym3_from_ym2_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       6, 7, 1)


## 4.1.25. Number of recaptures as young male 3 (7) from young male 1 (5) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_ym3_from_ym1_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       5, 7, 2)


## 4.1.26. Total number of recaptures as young male 4 (8) ----
# -------------------------------------------------------

total_nb_recaptures_ym4_true = nb.recaptures(lions.ch, 8)


## 4.1.27. Number of recaptures as young male 4 (8) at t + 1 ----
# ----------------------------------------------------------

total_nb_recaptures_ym4_t1_true = nb.recaptures.tn(lions.ch, 8, 1)


## 4.1.28. Number of recaptures as young male 4 (8) at t + 2 ----
# ----------------------------------------------------------

total_nb_recaptures_ym4_t2_true = nb.recaptures.tn(lions.ch, 8, 2)


## 4.1.29. Number of recaptures as young male 4 (8) from young male 3 (7) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_ym4_from_ym3_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                       7, 8, 2)


## 4.1.30. Number of recaptures as young male 4 (8) from young male 2 (6) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_ym4_from_ym2_t2_true = nb.recaptures.from.state.tn(lions.ch,
                                                                       6, 8, 2)

## 4.1.31. Total number of recaptures as nomads (9) ----
# -------------------------------------------------

total_nb_recaptures_nm_true = nb.recaptures(lions.ch, 9)


## 4.1.32. Number of recaptures as nomads (9) at t + 1 ----
# ----------------------------------------------------

total_nb_recaptures_nm_t1_true = nb.recaptures.tn(lions.ch, 9, 1)


## 4.1.33. Number of recaptures as nomads (9) at t + 2 ----
# ----------------------------------------------------

total_nb_recaptures_nm_t2_true = nb.recaptures.tn(lions.ch, 9, 2)


## 4.1.34. Number of recaptures as nomads (9) from male old subadult (3) at t + 2 ----
# -------------------------------------------------------------------------------

total_nb_recaptures_nm_from_sa2m_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       3, 9, 2)


## 4.1.35. Number of recaptures as nomads (9) from young male 1 (5) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym1_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                      5, 9, 1)


## 4.1.36. Number of recaptures as nomads (9) from young male 1 (5) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym1_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      5, 9, 2)


## 4.1.37. Number of recaptures as nomads (9) from young male 2 (6) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym2_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      6, 9, 1)


## 4.1.38. Number of recaptures as nomads (9) from young male 2 (6) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym2_t2_true = nb.recaptures.from.state.tn(lions.ch,
                                                                      6, 9, 2)


## 4.1.39. Number of recaptures as nomads (9) from young male 3 (7) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym3_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                      7, 9, 1)


## 4.1.40. Number of recaptures as nomads (9) from young male 3 (7) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym3_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      7, 9, 2)


## 4.1.41. Number of recaptures as nomads (9) from young male 4 (8) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym4_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      8, 9, 1)


## 4.1.42. Number of recaptures as nomads (9) from young male 4 (8) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym4_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      8, 9, 2)


## 4.1.43. Number of recaptures as nomads (9) from nomadic male (9) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_nm_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                     9, 9, 1)


## 4.1.44. Number of recaptures as nomads (9) from nomadic male (9) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_nm_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                     9, 9, 2)


## 4.1.45. Number of recaptures as nomads (9) from resident male (10) at t + 1 ----
# ----------------------------------------------------------------------------

total_nb_recaptures_nm_from_rm_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                     10, 9, 1)


## 4.1.46. Number of recaptures as nomads (9) from resident male (10) at t + 2 ----
# ----------------------------------------------------------------------------

total_nb_recaptures_nm_from_rm_t2_true = nb.recaptures.from.state.tn(lions.ch,
                                                                     10, 9, 2)


## 4.1.47. Total number of recaptures as residents (10) ----
# -----------------------------------------------------

total_nb_recaptures_rm_true = nb.recaptures(lions.ch, 10)


## 4.1.48. Number of recaptures as residents (10) at t + 1 ----
# --------------------------------------------------------

total_nb_recaptures_rm_t1_true = nb.recaptures.tn(lions.ch, 10, 1)


## 4.1.49. Number of recaptures as residents (10) at t + 2 ----
# --------------------------------------------------------

total_nb_recaptures_rm_t2_true = nb.recaptures.tn(lions.ch, 10, 2)


## 4.1.50. Number of recaptures as residents (10) from male old subadult (3) at t + 2 ----
# -----------------------------------------------------------------------------------

total_nb_recaptures_rm_from_sa2m_t2_true = nb.recaptures.from.state.tn(lions.ch,
                                                                       3, 10, 2)


## 4.1.51. Number of recaptures as residents (10) from young male 1 (5) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym1_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      5, 10, 1)


## 4.1.52. Number of recaptures as residents (10) from young male 1 (5) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym1_t2_true = nb.recaptures.from.state.tn(lions.ch,
                                                                      5, 10, 2)


## 4.1.53. Number of recaptures as residents (10) from young male 2 (6) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym2_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      6, 10, 1)


## 4.1.54. Number of recaptures as residents (10) from young male 2 (6) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym2_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      6, 10, 2)


## 4.1.55. Number of recaptures as residents (10) from young male 3 (7) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym3_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                      7, 10, 1)


## 4.1.56. Number of recaptures as residents (10) from young male 3 (7) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym3_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      7, 10, 2)


## 4.1.57. Number of recaptures as residents (10) from young male 4 (8) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym4_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      8, 10, 1)


## 4.1.58. Number of recaptures as residents (10) from young male 4 (8) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym4_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      8, 10, 2)


## 4.1.59. Number of recaptures as residents (10) from nomadic male (9) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_nm_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                     9, 10, 1)


## 4.1.60. Number of recaptures as residents (10) from nomadic male (9) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_nm_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                     9, 10, 2)


## 4.1.61. Number of recaptures as residents (10) from resident male (10) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_rm_t1_true = nb.recaptures.from.state.tn(lions.ch,
                                                                     10, 10, 1)


## 4.1.62. Number of recaptures as residents (10) from resident male (10) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_rm_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                     10, 10, 2)


## 4.1.63. Total number of recaptures as adult female (4) ----
# -------------------------------------------------------

total_nb_recaptures_af_true = nb.recaptures(lions.ch, 4)


## 4.1.64. Number of recaptures as adult female (4) at t + 1 ----
# ----------------------------------------------------------

total_nb_recaptures_af_t1_true = nb.recaptures.tn(lions.ch, 4, 1)


## 4.1.65. Number of recaptures as adult female (4) from female old subadults (2) at t + 1 ----
# ----------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa2f_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       2, 4, 1)


## 4.1.66. Number of recaptures as adult female (4) at t + 2 ----
# ----------------------------------------------------------

total_nb_recaptures_af_t2_true = nb.recaptures.tn(lions.ch, 4, 2)


## 4.1.67. Number of recaptures as adult female (4) from young subadults (1) at t + 2 ----
# -----------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa1_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                      1, 4, 2)


## 4.1.68. Number of recaptures as adult female (4) from female old subadults (2) at t + 2 ----
# ----------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa2f_t2_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                       2, 4, 2)


## 4.1.69. Number of recaptures as adult female (4) from adult female (4) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_af_from_af_t1_true = nb.recaptures.from.state.tn(lions.ch, 
                                                                     4, 4, 1)


## 4.1.70. Number of recaptures as adult female (4) from adult female (4) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_af_from_af_t2_true = nb.recaptures.from.state.tn(lions.ch,
                                                                     4, 4, 2)


## 4.1.71. Total number of dead recoveries (11) ----
# ---------------------------------------------

total_nb_dead_recoveries_true = nb.recaptures(lions.ch, 11)


## 4.2. Simulated data ----
# --------------------

# We apply the functions to each simulated dataset.

## 4.2.1. Total number of recaptures ----
# ----------------------------------

total_nb_recaptures_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_simulated = c(total_nb_recaptures_simulated, 
                                    nb.recaptures(data_simulation_output[[nsim]],
                                                  seq(1, 10)))
  
}  


## 4.2.2. Total number of recaptures at t + 1 ----
# -------------------------------------------

total_nb_recaptures_t1_simulated = vector(mode = "list",
                                          length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]], 
                                                              seq(1, 10), 1)
  
}


## 4.2.3. Total number of recaptures at t + 2 ----
# -------------------------------------------

total_nb_recaptures_t2_simulated = vector(mode = "list",
                                          length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                              seq(1, 10), 2)
  
}


## 4.2.4. Total number of recaptures as female old subadults (2) ----
# --------------------------------------------------------------

total_nb_recaptures_sa2f_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_sa2f_simulated = c(total_nb_recaptures_sa2f_simulated,
                                         nb.recaptures(data_simulation_output[[nsim]],
                                                       2))
  
} 


## 4.2.5. Number of recaptures as female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------

total_nb_recaptures_sa2f_t1_simulated = vector(mode = "list",
                                               length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_sa2f_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]], 
                                                                   2, 1)
  
} 


## 4.2.6. Total number of recaptures as male old subadults (3) ----
# ------------------------------------------------------------

total_nb_recaptures_sa2m_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_sa2m_simulated = c(total_nb_recaptures_sa2m_simulated, 
                                         nb.recaptures(data_simulation_output[[nsim]], 
                                                       3))
  
} 


## 4.2.7. Number of recaptures as male old subadults (3) at t + 1 ----
# ---------------------------------------------------------------

total_nb_recaptures_sa2m_t1_simulated = vector(mode = "list",
                                               length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_sa2m_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]], 
                                                                   3, 1)
  
}


## 4.2.8. Total number of recaptures as young males (5, 6, 7, 8) ----
# --------------------------------------------------------------

total_nb_recaptures_ym_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym_simulated = c(total_nb_recaptures_ym_simulated,
                                       nb.recaptures(data_simulation_output[[nsim]],
                                                     seq(5, 8)))
  
} 


## 4.2.9. Number of recaptures as young males (5, 6, 7, 8) at t + 1 ----
# -----------------------------------------------------------------

total_nb_recaptures_ym_t1_simulated = vector(mode = "list", 
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]], 
                                                                 seq(5, 8), 1)
  
} 


## 4.2.10. Number of recaptures as young males (5, 6, 7, 8) at t + 2 ----
# ------------------------------------------------------------------

total_nb_recaptures_ym_t2_simulated = vector(mode = "list",
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                 seq(5, 8), 2)
  
}


## 4.2.11. Total number of recaptures as young males 1 (5) ----
# --------------------------------------------------------

total_nb_recaptures_ym1_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym1_simulated = c(total_nb_recaptures_ym1_simulated, 
                                        nb.recaptures(data_simulation_output[[nsim]],
                                                      5))
  
} 


## 4.2.12. Number of recaptures as young males 1 (5) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym1_t1_simulated = vector(mode = "list", 
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym1_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  5, 1)
  
} 


## 4.2.13. Number of recaptures as young males 1 (5) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym1_t2_simulated = vector(mode = "list",
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym1_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  5, 2)
  
}


## 4.2.14. Number of recaptures as young males 1 (5) from male old subadult (3) at t + 1 ----
# --------------------------------------------------------------------------------------

total_nb_recaptures_ym1_from_sa2m_t1_simulated = vector(mode = "list",
                                                        length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym1_from_sa2m_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                       3, 5, 1)
  
} 


## 4.2.15. Number of recaptures as young males 1 (5) from young subadult (1) at t + 2 ----
# -----------------------------------------------------------------------------------

total_nb_recaptures_ym1_from_sa1_t2_simulated = vector(mode = "list", 
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym1_from_sa1_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                      1, 5, 2)
  
}


## 4.2.16. Total number of recaptures as young males 2 (6) ----
# --------------------------------------------------------

total_nb_recaptures_ym2_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym2_simulated = c(total_nb_recaptures_ym1_simulated,
                                        nb.recaptures(data_simulation_output[[nsim]],
                                                      6))
  
} 


## 4.2.17. Number of recaptures as young males 2 (6) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym2_t1_simulated = vector(mode = "list", 
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym2_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  6, 1)
  
} 


## 4.2.18. Number of recaptures as young males 2 (6) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym2_t2_simulated = vector(mode = "list", 
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym2_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  6, 2)
  
}


## 4.2.19. Number of recaptures as young males 2 (6) from young males 1 (5) at t + 1 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_ym2_from_ym1_t1_simulated = vector(mode = "list", 
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym2_from_ym1_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      5, 6, 1)
  
} 


## 4.2.20. Number of recaptures as young males 2 (6) from male old subadult (3) at t + 2 ----
# --------------------------------------------------------------------------------------

total_nb_recaptures_ym2_from_sa2m_t2_simulated = vector(mode = "list",
                                                        length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym2_from_sa2m_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                       3, 6, 2)
  
}


## 4.2.21. Total number of recaptures as young males 3 (7) ----
# --------------------------------------------------------

total_nb_recaptures_ym3_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym3_simulated = c(total_nb_recaptures_ym1_simulated, 
                                        nb.recaptures(data_simulation_output[[nsim]], 
                                                      7))
  
} 


## 4.2.22. Number of recaptures as young males 3 (7) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym3_t1_simulated = vector(mode = "list",
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym3_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  7, 1)
  
} 


## 4.2.23. Number of recaptures as young males 3 (7) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym3_t2_simulated = vector(mode = "list",
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym3_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  7, 2)
  
}


## 4.2.24. Number of recaptures as young males 3 (7) from young males 2 (6) at t + 1 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_ym3_from_ym2_t1_simulated = vector(mode = "list",
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym3_from_ym2_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      6, 7, 1)
  
} 


## 4.2.25. Number of recaptures as young males 3 (7) from young male 1 (5) at t + 2 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym3_from_ym1_t2_simulated = vector(mode = "list", 
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym3_from_ym1_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      5, 7, 2)
  
}


## 4.2.26. Total number of recaptures as young males 4 (8) ----
# --------------------------------------------------------

total_nb_recaptures_ym4_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym4_simulated = c(total_nb_recaptures_ym4_simulated, 
                                        nb.recaptures(data_simulation_output[[nsim]],
                                                      8))
  
} 


## 4.2.27. Number of recaptures as young males 4 (8) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym4_t1_simulated = vector(mode = "list",
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym4_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  8, 1)
  
} 


## 4.2.28. Number of recaptures as young males 4 (8) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym4_t2_simulated = vector(mode = "list",
                                              length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym4_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                  8, 2)
  
}


## 4.2.29. Number of recaptures as young males 4 (8) from young male 3 (7) at t + 1 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym4_from_ym3_t1_simulated = vector(mode = "list",
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym4_from_ym3_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      7, 8, 1)
  
} 


## 4.2.30. Number of recaptures as young males 4 (8) from young male 2 (6) at t + 2 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym4_from_ym2_t2_simulated = vector(mode = "list",
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_ym4_from_ym2_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                      6, 8, 2)
  
}


## 4.2.31. Total number of recaptures as nomads (9) ----
# -------------------------------------------------

total_nb_recaptures_nm_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_simulated = c(total_nb_recaptures_nm_simulated,
                                       nb.recaptures(data_simulation_output[[nsim]],
                                                     9))
  
} 


## 4.2.32. Number of recaptures as nomads (9) at t + 1 ----
# ----------------------------------------------------

total_nb_recaptures_nm_t1_simulated = vector(mode = "list",
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                 9, 1)
  
} 


## 4.2.33. Number of recaptures as nomads (9) at t + 2 ----
# ----------------------------------------------------

total_nb_recaptures_nm_t2_simulated = vector(mode = "list",
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                 9, 2)
  
}


## 4.2.34. Number of recaptures as nomads (9) from male old subadult (3) at t + 2 ----
# -------------------------------------------------------------------------------

total_nb_recaptures_nm_from_sa2m_t2_simulated = vector(mode = "list", 
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_sa2m_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      3, 9, 2)
  
}


## 4.2.35. Number of recaptures as nomads (9) from young male 1 (5) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym1_t1_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym1_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     5, 9, 1)
  
} 


## 4.2.36. Number of recaptures as nomads (9) from young male 1 (5) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym1_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym1_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                     5, 9, 2)
  
}


## 4.2.37. Number of recaptures as nomads (9) from young male 2 (6) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym2_t1_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym2_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                     6, 9, 1)
  
} 


## 4.2.38. Number of recaptures as nomads (9) from young male 2 (6) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym2_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym2_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     6, 9, 2)
  
}


## 4.2.39. Number of recaptures as nomads (9) from young male 3 (7) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym3_t1_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym3_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                     7, 9, 1)
  
} 


## 4.2.40. Number of recaptures as nomads (9) from young male 3 (7) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym3_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym3_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                     7, 9, 2)
  
}


## 4.2.41. Number of recaptures as nomads (9) from young male 4 (8) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym4_t1_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym4_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                     8, 9, 1)
  
} 


## 4.2.42. Number of recaptures as nomads (9) from young male 4 (8) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym4_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_ym4_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     8, 9, 2)
  
}


## 4.2.43. Number of recaptures as nomads (9) from nomadic male (9) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_nm_t1_simulated = vector(mode = "list",
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_nm_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                    9, 9, 1)
  
} 


## 4.2.44. Number of recaptures as nomads (9) from nomadic male (9) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_nm_t2_simulated = vector(mode = "list",
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_nm_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                    9, 9, 2)
  
}


## 4.2.45. Number of recaptures as nomads (9) from resident male (10) at t + 1 ----
# ----------------------------------------------------------------------------

total_nb_recaptures_nm_from_rm_t1_simulated = vector(mode = "list",
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_rm_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                    10, 9, 1)
  
} 


## 4.2.46. Number of recaptures as nomads (9) from resident male (10) at t + 2 ----
# ----------------------------------------------------------------------------

total_nb_recaptures_nm_from_rm_t2_simulated = vector(mode = "list",
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_nm_from_rm_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                    10, 9, 2)
  
}


## 4.2.47. Total number of recaptures as residents (10) ----
# -----------------------------------------------------

total_nb_recaptures_rm_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_simulated = c(total_nb_recaptures_rm_simulated, 
                                       nb.recaptures(data_simulation_output[[nsim]], 
                                                     10))
  
} 


## 4.2.48. Number of recaptures as residents (10) at t + 1 ----
# --------------------------------------------------------

total_nb_recaptures_rm_t1_simulated = vector(mode = "list", 
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]], 
                                                                 10, 1)
  
} 


## 4.2.49. Number of recaptures as residents (10) at t + 2 ----
# --------------------------------------------------------

total_nb_recaptures_rm_t2_simulated = vector(mode = "list", 
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                 10, 2)
  
}


## 4.2.50. Number of recaptures as residents (10) from male old subadult (3) at t + 2 ----
# -----------------------------------------------------------------------------------

total_nb_recaptures_rm_from_sa2m_t2_simulated = vector(mode = "list", 
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_sa2m_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                      3, 10, 2)
  
}


## 4.2.51. Number of recaptures as residents (10) from young male 1 (5) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym1_t1_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym1_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     5, 10, 1)
  
} 


## 4.2.52. Number of recaptures as residents (10) from young male 1 (5) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym1_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym1_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     5, 10, 2)
  
}


## 4.2.53. Number of recaptures as residents (10) from young male 2 (6) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym2_t1_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym2_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     6, 10, 1)
  
} 


## 4.2.54. Number of recaptures as residents (10) from young male 2 (6) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym2_t2_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym2_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     6, 10, 2)
  
}


## 4.2.55. Number of recaptures as residents (10) from young male 3 (7) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym3_t1_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym3_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     7, 10, 1)
  
} 


## 4.2.56. Number of recaptures as residents (10) from young male 3 (7) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym3_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym3_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     7, 10, 2)
  
}


## 4.2.57. Number of recaptures as residents (10) from young male 4 (8) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym4_t1_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym4_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     8, 10, 1)
  
} 


## 4.2.58. Number of recaptures as residents (10) from young male 4 (8) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym4_t2_simulated = vector(mode = "list", 
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_ym4_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     8, 10, 2)
  
}


## 4.2.59. Number of recaptures as residents (10) from nomadic male (9) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_nm_t1_simulated = vector(mode = "list", 
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_nm_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                    9, 10, 1)
  
} 


## 4.2.60. Number of recaptures as residents (10) from nomadic male (9) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_nm_t2_simulated = vector(mode = "list",
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_nm_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                    9, 10, 2)
  
}


## 4.2.61. Number of recaptures as residents (10) from resident male (10) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_rm_t1_simulated = vector(mode = "list", 
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_rm_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                    10, 10, 1)
  
} 


## 4.2.62. Number of recaptures as residents (10) from resident male (10) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_rm_t2_simulated = vector(mode = "list", 
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_rm_from_rm_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]], 
                                                                                    10, 10, 2)
  
}


## 4.2.63. Total number of recaptures as adult females (4) ----
# --------------------------------------------------------

total_nb_recaptures_af_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_simulated = c(total_nb_recaptures_af_simulated,
                                       nb.recaptures(data_simulation_output[[nsim]],
                                                     4))
  
} 


## 4.2.64. Number of recaptures as adult females (4) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_af_t1_simulated = vector(mode = "list", 
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_t1_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                 4, 1)
  
}


## 4.2.65. Number of recaptures as adult females (4) from female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa2f_t1_simulated = vector(mode = "list",
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_from_sa2f_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      2, 4, 1)
  
}


## 4.2.66. Number of recaptures as adult females (4) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_af_t2_simulated = vector(mode = "list", 
                                             length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_t2_simulated[[nsim]] = nb.recaptures.tn(data_simulation_output[[nsim]],
                                                                 4, 2)
  
}


## 4.2.67. Number of recaptures as adult females (4) from young subadults (1) at t + 2 ----
# ------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa1_t2_simulated = vector(mode = "list",
                                                      length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_from_sa1_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                     1, 4, 2)
  
}


## 4.2.68. Number of recaptures as adult females (4) from female old subadults (2) at t + 2 ----
# -----------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa2f_t2_simulated = vector(mode = "list",
                                                       length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_from_sa2f_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                      2, 4, 2)
  
}


## 4.2.69. Number of recaptures as adult females (4) from adult females (4) at t + 1 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_af_from_af_t1_simulated = vector(mode = "list",
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_from_af_t1_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                    4, 4, 1)
  
}


## 4.2.70. Number of recaptures as adult females (4) from adult females (4) at t + 2 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_af_from_af_t2_simulated = vector(mode = "list", 
                                                     length = length(data_simulation_output))

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_recaptures_af_from_af_t2_simulated[[nsim]] = nb.recaptures.from.state.tn(data_simulation_output[[nsim]],
                                                                                    4, 4, 2)
  
}


## 4.2.71. Total number of dead recoveries (11) ----
# ---------------------------------------------

total_nb_dead_recoveries_simulated = c()

for(nsim in 1:length(data_simulation_output)){
  
  total_nb_dead_recoveries_simulated = c(total_nb_dead_recoveries_simulated, 
                                         nb.recaptures(data_simulation_output[[nsim]],
                                                       11))
  
}




###########################################################################
#
# 5. Calculating Bayesian p-values ----
#
###########################################################################

## 5.1. Bayesian p-values ----
# -----------------------

# We calculate the Bayesian p-values as the proportion of simulate capture
# histories in which the metric value is above the observed value.

## 5.1.1. Total number of recaptures ----
# ----------------------------------

total_nb_recaptures_pvalue = length(which(total_nb_recaptures_simulated > total_nb_recaptures_true)) / length(data_simulation_output)


## 5.1.2. Total number of recaptures at t + 1 ----
# -------------------------------------------

total_nb_recaptures_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_t1_pvalue = c(total_nb_recaptures_t1_pvalue, 
                                    length(which(unlist(lapply(total_nb_recaptures_t1_simulated,
                                                               FUN = function(x) x[t] > total_nb_recaptures_t1_true[t])))) 
                                    / length(data_simulation_output))
}


## 5.1.3. Total number of recaptures at t + 2 ----
# -------------------------------------------

total_nb_recaptures_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_t2_pvalue = c(total_nb_recaptures_t2_pvalue,
                                    length(which(unlist(lapply(total_nb_recaptures_t2_simulated, 
                                                               FUN = function(x) x[t] > total_nb_recaptures_t2_true[t])))) 
                                    / length(data_simulation_output))
}


## 5.1.4. Total number of recaptures as female old subadults (2) ----
# --------------------------------------------------------------

total_nb_recaptures_sa2f_pvalue = length(which(total_nb_recaptures_sa2f_simulated > total_nb_recaptures_sa2f_true)) / length(data_simulation_output)


## 5.1.5. Number of recaptures as female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------

total_nb_recaptures_sa2f_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_sa2f_t1_pvalue = c(total_nb_recaptures_sa2f_t1_pvalue,
                                         length(which(unlist(lapply(total_nb_recaptures_sa2f_t1_simulated,
                                                                    FUN = function(x) x[t] > total_nb_recaptures_sa2f_t1_true[t])))) 
                                         / length(data_simulation_output))
}


## 5.1.6. Total number of recaptures as male old subadults (3) ----
# ------------------------------------------------------------

total_nb_recaptures_sa2m_pvalue = length(which(total_nb_recaptures_sa2m_simulated > total_nb_recaptures_sa2m_true)) / length(data_simulation_output)


## 5.1.7. Number of recaptures as male old subadults (3) at t + 1 ----
# ---------------------------------------------------------------

total_nb_recaptures_sa2m_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_sa2m_t1_pvalue = c(total_nb_recaptures_sa2m_t1_pvalue, 
                                         length(which(unlist(lapply(total_nb_recaptures_sa2m_t1_simulated, 
                                                                    FUN = function(x) x[t] > total_nb_recaptures_sa2m_t1_true[t])))) 
                                         / length(data_simulation_output))
}


## 5.1.8. Total number of recaptures as young males (5, 6, 7, 8) ----
# --------------------------------------------------------------

total_nb_recaptures_ym_pvalue = length(which(total_nb_recaptures_ym_simulated > total_nb_recaptures_ym_true)) / length(data_simulation_output)


## 5.1.9. Number of recaptures as young males (5, 6, 7, 8) at t + 1 ----
# -----------------------------------------------------------------

total_nb_recaptures_ym_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym_t1_pvalue = c(total_nb_recaptures_ym_t1_pvalue,
                                       length(which(unlist(lapply(total_nb_recaptures_ym_t1_simulated, 
                                                                  FUN = function(x) x[t] > total_nb_recaptures_ym_t1_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.10. Number of recaptures as young males (5, 6, 7, 8) at t + 2 ----
# ------------------------------------------------------------------

total_nb_recaptures_ym_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym_t2_pvalue = c(total_nb_recaptures_ym_t2_pvalue,
                                       length(which(unlist(lapply(total_nb_recaptures_ym_t2_simulated, 
                                                                  FUN = function(x) x[t] > total_nb_recaptures_ym_t2_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.11. Total number of recaptures as young males 1 (5) ----
# --------------------------------------------------------

total_nb_recaptures_ym1_pvalue = length(which(total_nb_recaptures_ym1_simulated > total_nb_recaptures_ym1_true)) / length(data_simulation_output)


## 5.1.12. Number of recaptures as young males 1 (5) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym1_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym1_t1_pvalue = c(total_nb_recaptures_ym1_t1_pvalue,
                                        length(which(unlist(lapply(total_nb_recaptures_ym1_t1_simulated,
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym1_t1_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.13. Number of recaptures as young males 1 (5) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym1_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym1_t2_pvalue = c(total_nb_recaptures_ym1_t2_pvalue,
                                        length(which(unlist(lapply(total_nb_recaptures_ym1_t2_simulated, 
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym1_t2_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.14. Number of recaptures as young males 1 (4) from male old subadult (3) at t + 1 ----
# --------------------------------------------------------------------------------------

total_nb_recaptures_ym1_from_sa2m_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym1_from_sa2m_t1_pvalue = c(total_nb_recaptures_ym1_from_sa2m_t1_pvalue, 
                                                  length(which(unlist(lapply(total_nb_recaptures_ym1_from_sa2m_t1_simulated, 
                                                                             FUN = function(x) x[t] > total_nb_recaptures_ym1_from_sa2m_t1_true[t])))) 
                                                  / length(data_simulation_output))
}


## 5.1.15. Number of recaptures as young males 1 (5) from young subadult (1) at t + 2 ----
# -----------------------------------------------------------------------------------

total_nb_recaptures_ym1_from_sa1_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym1_from_sa1_t2_pvalue = c(total_nb_recaptures_ym1_from_sa1_t2_pvalue,
                                                 length(which(unlist(lapply(total_nb_recaptures_ym1_from_sa1_t2_simulated, 
                                                                            FUN = function(x) x[t] > total_nb_recaptures_ym1_from_sa1_t2_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.16. Total number of recaptures as young males 2 (6) ----
# --------------------------------------------------------

total_nb_recaptures_ym2_pvalue = length(which(total_nb_recaptures_ym2_simulated > total_nb_recaptures_ym2_true)) / length(data_simulation_output)


## 5.1.17. Number of recaptures as young males 2 (6) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym2_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym2_t1_pvalue = c(total_nb_recaptures_ym2_t1_pvalue,
                                        length(which(unlist(lapply(total_nb_recaptures_ym2_t1_simulated,
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym2_t1_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.18. Number of recaptures as young males 2 (6) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym2_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym2_t2_pvalue = c(total_nb_recaptures_ym2_t2_pvalue,
                                        length(which(unlist(lapply(total_nb_recaptures_ym2_t2_simulated, 
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym2_t2_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.19. Number of recaptures as young males 2 (6) from young male 1 (5) at t + 1 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym2_from_ym1_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym2_from_ym1_t1_pvalue = c(total_nb_recaptures_ym2_from_ym1_t1_pvalue, 
                                                 length(which(unlist(lapply(total_nb_recaptures_ym2_from_ym1_t1_simulated, 
                                                                            FUN = function(x) x[t] > total_nb_recaptures_ym2_from_ym1_t1_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.20. Number of recaptures as young males 2 (6) from male old subadult (3) at t + 2 ----
# --------------------------------------------------------------------------------------

total_nb_recaptures_ym2_from_sa2m_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym2_from_sa2m_t2_pvalue = c(total_nb_recaptures_ym2_from_sa2m_t2_pvalue, 
                                                  length(which(unlist(lapply(total_nb_recaptures_ym2_from_sa2m_t2_simulated, 
                                                                             FUN = function(x) x[t] > total_nb_recaptures_ym2_from_sa2m_t2_true[t])))) 
                                                  / length(data_simulation_output))
}


## 5.1.21. Total number of recaptures as young males 3 (7) ----
# --------------------------------------------------------

total_nb_recaptures_ym3_pvalue = length(which(total_nb_recaptures_ym3_simulated > total_nb_recaptures_ym3_true)) / length(data_simulation_output)


## 5.1.22. Number of recaptures as young males 3 (7) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym3_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym3_t1_pvalue = c(total_nb_recaptures_ym3_t1_pvalue, 
                                        length(which(unlist(lapply(total_nb_recaptures_ym3_t1_simulated, 
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym3_t1_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.23. Number of recaptures as young males 3 (7) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym3_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym3_t2_pvalue = c(total_nb_recaptures_ym3_t2_pvalue,
                                        length(which(unlist(lapply(total_nb_recaptures_ym3_t2_simulated, 
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym3_t2_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.24. Number of recaptures as young males 3 (7) from young male 2 (6) at t + 1 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym3_from_ym2_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym3_from_ym2_t1_pvalue = c(total_nb_recaptures_ym3_from_ym2_t1_pvalue,
                                                 length(which(unlist(lapply(total_nb_recaptures_ym3_from_ym2_t1_simulated,
                                                                            FUN = function(x) x[t] > total_nb_recaptures_ym3_from_ym2_t1_true[t]))))
                                                 / length(data_simulation_output))
}


## 5.1.25. Number of recaptures as young males 3 (7) from young male 1 (5) at t + 2 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym3_from_ym1_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym3_from_ym1_t2_pvalue = c(total_nb_recaptures_ym3_from_ym1_t2_pvalue, 
                                                 length(which(unlist(lapply(total_nb_recaptures_ym3_from_ym1_t2_simulated, 
                                                                            FUN = function(x) x[t] > total_nb_recaptures_ym3_from_ym1_t2_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.26. Total number of recaptures as young males 4 (8) ----
# --------------------------------------------------------

total_nb_recaptures_ym4_pvalue = length(which(total_nb_recaptures_ym4_simulated > total_nb_recaptures_ym4_true)) / length(data_simulation_output)


## 5.1.27. Number of recaptures as young males 4 (8) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_ym4_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym4_t1_pvalue = c(total_nb_recaptures_ym4_t1_pvalue, 
                                        length(which(unlist(lapply(total_nb_recaptures_ym4_t1_simulated,
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym4_t1_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.28. Number of recaptures as young males 4 (8) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_ym4_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym4_t2_pvalue = c(total_nb_recaptures_ym4_t2_pvalue,
                                        length(which(unlist(lapply(total_nb_recaptures_ym4_t2_simulated,
                                                                   FUN = function(x) x[t] > total_nb_recaptures_ym4_t2_true[t])))) 
                                        / length(data_simulation_output))
}


## 5.1.29. Number of recaptures as young males 4 (8) from young male 3 (7) at t + 1 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym4_from_ym3_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_ym4_from_ym3_t1_pvalue = c(total_nb_recaptures_ym4_from_ym3_t1_pvalue, 
                                                 length(which(unlist(lapply(total_nb_recaptures_ym4_from_ym3_t1_simulated,
                                                                            FUN = function(x) x[t] > total_nb_recaptures_ym4_from_ym3_t1_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.30. Number of recaptures as young males 4 (8) from young male 2 (6) at t + 2 ----
# ---------------------------------------------------------------------------------

total_nb_recaptures_ym4_from_ym2_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_ym4_from_ym2_t2_pvalue = c(total_nb_recaptures_ym4_from_ym2_t2_pvalue,
                                                 length(which(unlist(lapply(total_nb_recaptures_ym4_from_ym2_t2_simulated,
                                                                            FUN = function(x) x[t] > total_nb_recaptures_ym4_from_ym2_t2_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.31. Total number of recaptures as nomads (9) ----
# -------------------------------------------------

total_nb_recaptures_nm_pvalue = length(which(total_nb_recaptures_nm_simulated > total_nb_recaptures_nm_true)) / length(data_simulation_output)


## 5.1.32. Number of recaptures as nomads (9) at t + 1 ----
# ----------------------------------------------------

total_nb_recaptures_nm_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_t1_pvalue = c(total_nb_recaptures_nm_t1_pvalue,
                                       length(which(unlist(lapply(total_nb_recaptures_nm_t1_simulated, 
                                                                  FUN = function(x) x[t] > total_nb_recaptures_nm_t1_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.33. Number of recaptures as nomads (9) at t + 2 ----
# ----------------------------------------------------

total_nb_recaptures_nm_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_t2_pvalue = c(total_nb_recaptures_nm_t2_pvalue,
                                       length(which(unlist(lapply(total_nb_recaptures_nm_t2_simulated,
                                                                  FUN = function(x) x[t] > total_nb_recaptures_nm_t2_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.34. Number of recaptures as nomads (9) from male old subadult (3) at t + 2 ----
# -------------------------------------------------------------------------------

total_nb_recaptures_nm_from_sa2m_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_sa2m_t2_pvalue = c(total_nb_recaptures_nm_from_sa2m_t2_pvalue,
                                                 length(which(unlist(lapply(total_nb_recaptures_nm_from_sa2m_t2_simulated,
                                                                            FUN = function(x) x[t] > total_nb_recaptures_nm_from_sa2m_t2_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.35. Number of recaptures as nomads (9) from young male 1 (5) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym1_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_from_ym1_t1_pvalue = c(total_nb_recaptures_nm_from_ym1_t1_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym1_t1_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym1_t1_true[t]))))
                                                / length(data_simulation_output))
}


## 5.1.36. Number of recaptures as nomads (9) from young male 1 (5) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym1_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_ym1_t2_pvalue = c(total_nb_recaptures_nm_from_ym1_t2_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym1_t2_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym1_t2_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.37. Number of recaptures as nomads (9) from young male 2 (6) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym2_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_from_ym2_t1_pvalue = c(total_nb_recaptures_nm_from_ym2_t1_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym2_t1_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym2_t1_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.38. Number of recaptures as nomads (9) from young male 2 (6) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym2_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_ym2_t2_pvalue = c(total_nb_recaptures_nm_from_ym2_t2_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym2_t2_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym2_t2_true[t]))))
                                                / length(data_simulation_output))
}


## 5.1.39. Number of recaptures as nomads (9) from young male 3 (7) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym3_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_from_ym3_t1_pvalue = c(total_nb_recaptures_nm_from_ym3_t1_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym3_t1_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym3_t1_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.40. Number of recaptures as nomads (9) from young male 3 (7) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym3_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_ym3_t2_pvalue = c(total_nb_recaptures_nm_from_ym3_t2_pvalue, 
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym3_t2_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym3_t2_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.41. Number of recaptures as nomads (9) from young male 4 (8) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym4_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_from_ym4_t1_pvalue = c(total_nb_recaptures_nm_from_ym4_t1_pvalue, 
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym4_t1_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym4_t1_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.42. Number of recaptures as nomads (9) from young male 4 (8) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_ym4_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_ym4_t2_pvalue = c(total_nb_recaptures_nm_from_ym4_t2_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_nm_from_ym4_t2_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_nm_from_ym4_t2_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.43. Number of recaptures as nomads (9) from nomadic male (9) at t + 1 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_nm_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_from_nm_t1_pvalue = c(total_nb_recaptures_nm_from_nm_t1_pvalue, 
                                               length(which(unlist(lapply(total_nb_recaptures_nm_from_nm_t1_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_nm_from_nm_t1_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.44. Number of recaptures as nomads (9) from nomadic male (9) at t + 2 ----
# --------------------------------------------------------------------------

total_nb_recaptures_nm_from_nm_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_nm_t2_pvalue = c(total_nb_recaptures_nm_from_nm_t2_pvalue, 
                                               length(which(unlist(lapply(total_nb_recaptures_nm_from_nm_t2_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_nm_from_nm_t2_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.45. Number of recaptures as nomads (9) from resident male (10) at t + 1 ----
# ----------------------------------------------------------------------------

total_nb_recaptures_nm_from_rm_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_nm_from_rm_t1_pvalue = c(total_nb_recaptures_nm_from_rm_t1_pvalue, 
                                               length(which(unlist(lapply(total_nb_recaptures_nm_from_rm_t1_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_nm_from_rm_t1_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.46. Number of recaptures as nomads (9) from resident male (10) at t + 2 ----
# ----------------------------------------------------------------------------

total_nb_recaptures_nm_from_rm_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_nm_from_rm_t2_pvalue = c(total_nb_recaptures_nm_from_rm_t2_pvalue,
                                               length(which(unlist(lapply(total_nb_recaptures_nm_from_rm_t2_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_nm_from_rm_t2_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.47. Total number of recaptures as residents (10) ----
# -----------------------------------------------------

total_nb_recaptures_rm_pvalue = length(which(total_nb_recaptures_rm_simulated > total_nb_recaptures_rm_true)) / length(data_simulation_output)


## 5.1.48. Number of recaptures as residents (10) at t + 1 ----
# -------------------------------------------------------

total_nb_recaptures_rm_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_t1_pvalue = c(total_nb_recaptures_rm_t1_pvalue,
                                       length(which(unlist(lapply(total_nb_recaptures_rm_t1_simulated, 
                                                                  FUN = function(x) x[t] > total_nb_recaptures_rm_t1_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.49. Number of recaptures as residents (10) at t + 2 ----
# -------------------------------------------------------

total_nb_recaptures_rm_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_t2_pvalue = c(total_nb_recaptures_rm_t2_pvalue, 
                                       length(which(unlist(lapply(total_nb_recaptures_rm_t2_simulated,
                                                                  FUN = function(x) x[t] > total_nb_recaptures_rm_t2_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.50. Number of recaptures as residents (10) from male old subadult (3) at t + 2 ----
# -----------------------------------------------------------------------------------

total_nb_recaptures_rm_from_sa2m_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_sa2m_t2_pvalue = c(total_nb_recaptures_rm_from_sa2m_t2_pvalue,
                                                 length(which(unlist(lapply(total_nb_recaptures_rm_from_sa2m_t2_simulated, 
                                                                            FUN = function(x) x[t] > total_nb_recaptures_rm_from_sa2m_t2_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.51. Number of recaptures as residents (10) from young male 1 (5) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym1_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_from_ym1_t1_pvalue = c(total_nb_recaptures_rm_from_ym1_t1_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym1_t1_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym1_t1_true[t]))))
                                                / length(data_simulation_output))
}


## 5.1.52. Number of recaptures as residents (10) from young male 1 (5) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym1_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_ym1_t2_pvalue = c(total_nb_recaptures_rm_from_ym1_t2_pvalue, 
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym1_t2_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym1_t2_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.53. Number of recaptures as residents (10) from young male 2 (6) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym2_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_from_ym2_t1_pvalue = c(total_nb_recaptures_rm_from_ym2_t1_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym2_t1_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym2_t1_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.54. Number of recaptures as residents (10) from young male 2 (6) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym2_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_ym2_t2_pvalue = c(total_nb_recaptures_rm_from_ym2_t2_pvalue, 
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym2_t2_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym2_t2_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.55. Number of recaptures as residents (10) from young male 3 (7) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym3_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_from_ym3_t1_pvalue = c(total_nb_recaptures_rm_from_ym3_t1_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym3_t1_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym3_t1_true[t]))))
                                                / length(data_simulation_output))
}


## 5.1.56. Number of recaptures as residents (10) from young male 3 (7) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym3_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_ym3_t2_pvalue = c(total_nb_recaptures_rm_from_ym3_t2_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym3_t2_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym3_t2_true[t]))))
                                                / length(data_simulation_output))
}


## 5.1.57. Number of recaptures as residents (10) from young male 4 (8) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym4_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_from_ym4_t1_pvalue = c(total_nb_recaptures_rm_from_ym4_t1_pvalue, 
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym4_t1_simulated,
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym4_t1_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.58. Number of recaptures as residents (10) from young male 4 (8) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_ym4_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_ym4_t2_pvalue = c(total_nb_recaptures_rm_from_ym4_t2_pvalue, 
                                                length(which(unlist(lapply(total_nb_recaptures_rm_from_ym4_t2_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_rm_from_ym4_t2_true[t]))))
                                                / length(data_simulation_output))
}


## 5.1.59. Number of recaptures as residents (10) from nomadic male (9) at t + 1 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_nm_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_from_nm_t1_pvalue = c(total_nb_recaptures_rm_from_nm_t1_pvalue, 
                                               length(which(unlist(lapply(total_nb_recaptures_rm_from_nm_t1_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_rm_from_nm_t1_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.60. Number of recaptures as residents (10) from nomadic male (9) at t + 2 ----
# ------------------------------------------------------------------------------

total_nb_recaptures_rm_from_nm_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_nm_t2_pvalue = c(total_nb_recaptures_rm_from_nm_t2_pvalue,
                                               length(which(unlist(lapply(total_nb_recaptures_rm_from_nm_t2_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_rm_from_nm_t2_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.61. Number of recaptures as residents (10) from resident male (10) at t + 1 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_rm_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_rm_from_rm_t1_pvalue = c(total_nb_recaptures_rm_from_rm_t1_pvalue, 
                                               length(which(unlist(lapply(total_nb_recaptures_rm_from_rm_t1_simulated, 
                                                                          FUN = function(x) x[t] > total_nb_recaptures_rm_from_rm_t1_true[t]))))
                                               / length(data_simulation_output))
}


## 5.1.62. Number of recaptures as residents (10) from resident male (10) at t + 2 ----
# --------------------------------------------------------------------------------

total_nb_recaptures_rm_from_rm_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_rm_from_rm_t2_pvalue = c(total_nb_recaptures_rm_from_rm_t2_pvalue,
                                               length(which(unlist(lapply(total_nb_recaptures_rm_from_rm_t2_simulated, 
                                                                          FUN = function(x) x[t] > total_nb_recaptures_rm_from_rm_t2_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.63. Total number of recaptures as adult females (4) ----
# --------------------------------------------------------

total_nb_recaptures_af_pvalue = length(which(total_nb_recaptures_af_simulated > total_nb_recaptures_af_true)) / length(data_simulation_output)


## 5.1.64. Number of recaptures as adult females (4) at t + 1 ----
# -----------------------------------------------------------

total_nb_recaptures_af_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_af_t1_pvalue = c(total_nb_recaptures_af_t1_pvalue, 
                                       length(which(unlist(lapply(total_nb_recaptures_af_t1_simulated,
                                                                  FUN = function(x) x[t] > total_nb_recaptures_af_t1_true[t])))) 
                                       / length(data_simulation_output))
}


## 5.1.65. Number of recaptures as adult females (4) from female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa2f_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_af_from_sa2f_t1_pvalue = c(total_nb_recaptures_af_from_sa2f_t1_pvalue, 
                                                 length(which(unlist(lapply(total_nb_recaptures_af_from_sa2f_t1_simulated, 
                                                                            FUN = function(x) x[t] > total_nb_recaptures_af_from_sa2f_t1_true[t]))))
                                                 / length(data_simulation_output))
}


## 5.1.66. Number of recaptures as adult females (4) at t + 2 ----
# -----------------------------------------------------------

total_nb_recaptures_af_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 2)){
  
  total_nb_recaptures_af_t2_pvalue = c(total_nb_recaptures_af_t2_pvalue,
                                       length(which(unlist(lapply(total_nb_recaptures_af_t2_simulated, 
                                                                  FUN = function(x) x[t] > total_nb_recaptures_af_t2_true[t]))))
                                       / length(data_simulation_output))
}


## 5.1.67. Number of recaptures as adult females (4) from young subadults (1) at t + 2 ----
# ------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa1_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_af_from_sa1_t2_pvalue = c(total_nb_recaptures_af_from_sa1_t2_pvalue,
                                                length(which(unlist(lapply(total_nb_recaptures_af_from_sa1_t2_simulated, 
                                                                           FUN = function(x) x[t] > total_nb_recaptures_af_from_sa1_t2_true[t])))) 
                                                / length(data_simulation_output))
}


## 5.1.68. Number of recaptures as adult females (4) from female old subadults (2) at t + 2 ----
# -----------------------------------------------------------------------------------------

total_nb_recaptures_af_from_sa2f_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_af_from_sa2f_t2_pvalue = c(total_nb_recaptures_af_from_sa2f_t2_pvalue, 
                                                 length(which(unlist(lapply(total_nb_recaptures_af_from_sa2f_t2_simulated,
                                                                            FUN = function(x) x[t] > total_nb_recaptures_af_from_sa2f_t2_true[t])))) 
                                                 / length(data_simulation_output))
}


## 5.1.69. Number of recaptures as adult females (4) from adult females (4) at t + 1 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_af_from_af_t1_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_af_from_af_t1_pvalue = c(total_nb_recaptures_af_from_af_t1_pvalue,
                                               length(which(unlist(lapply(total_nb_recaptures_af_from_af_t1_simulated, 
                                                                          FUN = function(x) x[t] > total_nb_recaptures_af_from_af_t1_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.70. Number of recaptures as adult females (4) from adult females (4) at t + 2 ----
# ----------------------------------------------------------------------------------

total_nb_recaptures_af_from_af_t2_pvalue = c()

for(t in 1:(ncol(lions.ch) - 1)){
  
  total_nb_recaptures_af_from_af_t2_pvalue = c(total_nb_recaptures_af_from_af_t2_pvalue,
                                               length(which(unlist(lapply(total_nb_recaptures_af_from_af_t2_simulated,
                                                                          FUN = function(x) x[t] > total_nb_recaptures_af_from_af_t2_true[t])))) 
                                               / length(data_simulation_output))
}


## 5.1.71. Total number of dead recoveries (11) ----
# ---------------------------------------------

total_nb_dead_recoveries_pvalue = length(which(total_nb_dead_recoveries_simulated > total_nb_dead_recoveries_true)) / length(data_simulation_output)




###########################################################################
#
# 6. Plotting simulated and true metric values and Bayesian p-values ----
#
###########################################################################

## 6.1. Plotting simulated and true values ----
# ----------------------------------------

# Plot settings:

# Define background color of plot 
colBG = "transparent"

# Define color, font, font size of plot
colPlot = "black"
font = "Helvetica"
fontSize = 7


# We first merge all the simulated and observed values in the same
# dataframe before plotting everything together.

# We do not plot the metrics for which the observed and simulated
# values are 0, that is:
#
# Number of recaptures as young male 4 (8) at t + 1 
# Number of recaptures as young male 4 (8) from young male 3 (7) at t + 1
# Number of recaptures as nomads (9) from young male 3 (7) at t + 1
# Number of recaptures as nomads (9) from young male 3 (7) at t + 2
# Number of recaptures as nomads (9) from young male 4 (8) at t + 1
# Number of recaptures as nomads (9) from young male 4 (8) at t + 2
# Number of recaptures as residents (10) from young male 3 (7) at t + 1
# Number of recaptures as residents (10) from young male 3 (7) at t + 2
# Number of recaptures as residents (10) from young male 4 (8) at t + 1
# Number of recaptures as residents (10) from young male 4 (8) at t + 2


## 6.1.1. Total number of recaptures ----
# ----------------------------------

plot_title = "Total nb recaptures"
sim_distribution = total_nb_recaptures_simulated
true_value = total_nb_recaptures_true

sim_obs_df = data.frame(metric = plot_title,
                        sim_distribution = sim_distribution,
                        true_value = true_value)


## 6.1.2. Total number of recaptures at t + 1 ----
# -------------------------------------------

plot_title = "Total nb recaptures\nat t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_t1_simulated, sum))
true_value = sum(total_nb_recaptures_t1_true)

sim_obs_df = rbind(sim_obs_df, 
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.3. Total number of recaptures at t + 2 ----
# -------------------------------------------

plot_title = "Total nb recaptures\nat t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_t2_simulated, sum))
true_value = sum(total_nb_recaptures_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.4. Total number of recaptures as female old subadults (2) ----
# --------------------------------------------------------------

plot_title = "Nb recaptures as\nfemale old subadults"
sim_distribution = total_nb_recaptures_sa2f_simulated
true_value = total_nb_recaptures_sa2f_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.5. Number of recaptures as female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------

plot_title = "Nb recaptures as\nfemale old subadults\nat t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_sa2f_t1_simulated, sum))
true_value = sum(total_nb_recaptures_sa2f_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.6. Total number of recaptures as male old subadults (3) ----
# ------------------------------------------------------------

plot_title = "Nb recaptures as\nmale old subadults"
sim_distribution = total_nb_recaptures_sa2m_simulated
true_value = total_nb_recaptures_sa2m_true

sim_obs_df = rbind(sim_obs_df, 
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.7. Number of recaptures as male old subadults (3) at t + 1 ----
# ---------------------------------------------------------------

plot_title = "Nb recaptures as male\nold subadults at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_sa2m_t1_simulated, sum))
true_value = sum(total_nb_recaptures_sa2m_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.8. Total number of recaptures as young males (5, 6, 7, 8) ----
# --------------------------------------------------------------

plot_title = "Nb recaptures\nas young male"
sim_distribution = total_nb_recaptures_ym_simulated
true_value = total_nb_recaptures_ym_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.9. Number of recaptures as young males (5, 6, 7, 8) at t + 1 ----
# -----------------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.10. Number of recaptures as young males (5, 6, 7, 8) at t + 2 ----
# ------------------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.11. Total number of recaptures as young males 1 (5) ----
# --------------------------------------------------------

plot_title = "Nb recaptures\nas young male 1"
sim_distribution = total_nb_recaptures_ym1_simulated
true_value = total_nb_recaptures_ym1_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.12. Number of recaptures as young males 1 (5) at t + 1 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 1 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym1_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym1_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.13. Number of recaptures as young males 1 (5) at t + 2 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 1 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym1_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym1_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.14. Number of recaptures as young males 1 (5) from male old subadult (3) at t + 1 ----
# --------------------------------------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 1 from\nmale old subadult at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym1_from_sa2m_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym1_from_sa2m_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.15. Number of recaptures as young males 1 (5) from young subadult (1) at t + 2 ----
# -----------------------------------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 1 from\nyoung subadult at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym1_from_sa1_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym1_from_sa1_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.16. Total number of recaptures as young males 2 (6) ----
# -------------------------------------------------------

plot_title = "Nb recaptures\nas young male 2"
sim_distribution = total_nb_recaptures_ym2_simulated
true_value = total_nb_recaptures_ym2_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.17. Number of recaptures as young males 2 (6) at t + 1 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 2 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym2_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym2_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.18. Number of recaptures as young males 2 (6) at t + 2 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 2 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym2_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym2_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.19. Number of recaptures as young males 2 (6) from young male 1 (5) at t + 1 ----
# ---------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas young male 2 from\nyoung male 1 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym2_from_ym1_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym2_from_ym1_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.20. Number of recaptures as young males 2 (6) from male old subadult (3) at t + 2 ----
# --------------------------------------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 2 from male\nold subadult at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym2_from_sa2m_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym2_from_sa2m_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.21. Total number of recaptures as young males 3 (7) ----
# -------------------------------------------------------

plot_title = "Nb recaptures\nas young male 3"
sim_distribution = total_nb_recaptures_ym3_simulated
true_value = total_nb_recaptures_ym3_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.22. Number of recaptures as young males 3 (7) at t + 1 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 3 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym3_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym3_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.23. Number of recaptures as young males 3 (7) at t + 2 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 3 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym3_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym3_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.24. Number of recaptures as young males 3 (7) from young male 2 (6) at t + 1 ----
# ---------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas young male 3 from\nyoung male 2 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym3_from_ym2_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym3_from_ym2_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.25. Number of recaptures as young males 3 (7) from young male 1 (5) at t + 2 ----
# ---------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas young male 3 from\nyoung male 1 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym3_from_ym1_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym3_from_ym1_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.26. Total number of recaptures as young males 4 (8) ----
# -------------------------------------------------------

plot_title = "Nb recaptures\nas young male 4"
sim_distribution = total_nb_recaptures_ym4_simulated
true_value = total_nb_recaptures_ym4_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.27. Number of recaptures as young males 4 (8) at t + 1 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 4 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym4_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym4_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.28. Number of recaptures as young males 4 (8) at t + 2 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nyoung male 4 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym4_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym4_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.29. Number of recaptures as young males 4 (8) from young male 3 (7) at t + 1 ----
# ---------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas young male 4 from\nyoung male 3 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_ym4_from_ym3_t1_simulated, sum))
true_value = sum(total_nb_recaptures_ym4_from_ym3_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.30. Number of recaptures as young males 4 (8) from young male 2 (6) at t + 2 ----
# ---------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas young male 4 from\nyoung male 2 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_ym4_from_ym2_t2_simulated, sum))
true_value = sum(total_nb_recaptures_ym4_from_ym2_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.31. Total number of recaptures as nomads (9) ----
# -------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male"
sim_distribution = total_nb_recaptures_nm_simulated
true_value = total_nb_recaptures_nm_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.32. Number of recaptures as nomads (9) at t + 1 ----
# ----------------------------------------------------

plot_title = "Nb recaptures as\nnomadic male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.33. Number of recaptures as nomads (9) at t + 2 ----
# ----------------------------------------------------

plot_title = "Nb recaptures as\nnomadic male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.34. Number of recaptures as nomads (9) from male old subadult (3) at t + 2 ----
# -------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nmale old subadult at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_sa2m_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_sa2m_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.35. Number of recaptures as nomads (9) from young male 1 (5) at t + 1 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 1 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym1_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym1_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.36. Number of recaptures as nomads (9) from young male 1 (5) at t + 2 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 1 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym1_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym1_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.37. Number of recaptures as nomads (9) from young male 2 (6) at t + 1 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 2 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym2_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym2_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.38. Number of recaptures as nomads (9) from young male 2 (6) at t + 2 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 2 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym2_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym2_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.39. Number of recaptures as nomads (9) from young male 3 (7) at t + 1 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 3 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym3_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym3_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.40. Number of recaptures as nomads (9) from young male 3 (7) at t + 2 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 3 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym3_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym3_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.41. Number of recaptures as nomads (9) from young male 4 (8) at t + 1 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 4 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym4_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym4_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.42. Number of recaptures as nomads (9) from young male 4 (8) at t + 2 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nyoung male 4 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_ym4_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_ym4_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.43. Number of recaptures as nomads (9) from nomadic male (9) at t + 1 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nnomadic male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_nm_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_nm_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.44. Number of recaptures as nomads (9) from nomadic male (9) at t + 2 ----
# --------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nnomadic male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_nm_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_nm_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.45. Number of recaptures as nomads (9) from resident male (10) at t + 1 ----
# ----------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nresident male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_rm_t1_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_rm_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.46. Number of recaptures as nomads (9) from resident male (10) at t + 2 ----
# ----------------------------------------------------------------------------

plot_title = "Nb recaptures\nas nomadic male from\nresident male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_nm_from_nm_t2_simulated, sum))
true_value = sum(total_nb_recaptures_nm_from_nm_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.47. Total number of recaptures as residents (10) ----
# -----------------------------------------------------

plot_title = "Nb recaptures\nas resident male"
sim_distribution = total_nb_recaptures_rm_simulated
true_value = total_nb_recaptures_rm_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.48. Number of recaptures as residents (10) at t + 1 ----
# --------------------------------------------------------

plot_title = "Nb recaptures as\nresident male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.49. Number of recaptures as residents (10) at t + 2 ----
# --------------------------------------------------------

plot_title = "Nb recaptures as\nresident male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.50. Number of recaptures as residents (10) from male old subadult (3) at t + 2 ----
# -----------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nmale old subadult at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_sa2m_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_sa2m_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.51. Number of recaptures as residents (10) from young male 1 (4) at t + 1 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 1 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym1_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym1_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.52. Number of recaptures as residents (10) from young male 1 (5) at t + 2 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 1 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym1_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym1_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.53. Number of recaptures as residents (10) from young male 2 (6) at t + 1 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 2 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym2_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym2_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.54. Number of recaptures as residents (10) from young male 2 (6) at t + 2 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 2 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym2_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym2_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.55. Number of recaptures as residents (10) from young male 3 (7) at t + 1 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 3 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym3_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym3_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.56. Number of recaptures as residents (10) from young male 3 (7) at t + 2 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 3 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym3_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym3_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.57. Number of recaptures as residents (10) from young male 4 (8) at t + 1 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 4 at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym4_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym4_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.58. Number of recaptures as residents (10) from young male 4 (8) at t + 2 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nyoung male 4 at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_ym4_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_ym4_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.59. Number of recaptures as residents (10) from nomadic male (9) at t + 1 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nnomadic male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_nm_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_nm_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.60. Number of recaptures as residents (10) from nomadic male (9) at t + 2 ----
# ------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nnomadic male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_nm_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_nm_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.61. Number of recaptures as residents (10) from resident male (10) at t + 1 ----
# --------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nresident male at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_rm_t1_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_rm_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.62. Number of recaptures as residents (10) from resident male (10) at t + 2 ----
# --------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas resident male from\nresident male at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_rm_from_rm_t2_simulated, sum))
true_value = sum(total_nb_recaptures_rm_from_rm_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.63. Total number of recaptures as adult females (4) ----
# --------------------------------------------------------

plot_title = "Nb recaptures\nas adult female"
sim_distribution = total_nb_recaptures_af_simulated
true_value = total_nb_recaptures_af_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.64. Number of recaptures as adult females (4) at t + 1 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nadult female at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_af_t1_simulated, sum))
true_value = sum(total_nb_recaptures_af_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.65. Number of recaptures as adult females (4) from female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------------------------------

plot_title = "Nb recaptures as\nadult female from\nfemale old subadult\nat t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_af_from_sa2f_t1_simulated, sum))
true_value = sum(total_nb_recaptures_af_from_sa2f_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.66. Number of recaptures as adult females (4) at t + 2 ----
# -----------------------------------------------------------

plot_title = "Nb recaptures as\nadult female at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_af_t2_simulated, sum))
true_value = sum(total_nb_recaptures_af_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.67. Number of recaptures as adult females (4) from young subadults (1) at t + 2 ----
# ------------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas adult female from\nyoung subadult at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_af_from_sa1_t2_simulated, sum))
true_value = sum(total_nb_recaptures_af_from_sa1_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.68. Number of recaptures as adult females (4) from female old subadults (2) at t + 2 ----
# -----------------------------------------------------------------------------------------

plot_title = "Nb recaptures as\nadult female from\nfemale old subadult\nat t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_af_from_sa2f_t2_simulated, sum))
true_value = sum(total_nb_recaptures_af_from_sa2f_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.69. Number of recaptures as adult females (4) from adult females (4) at t + 1 ----
# ----------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas adult female from\nadult female at t+1"
sim_distribution = unlist(lapply(total_nb_recaptures_af_from_af_t1_simulated, sum))
true_value = sum(total_nb_recaptures_af_from_af_t1_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.70. Number of recaptures as adult females (4) from adult females (4) at t + 2 ----
# ----------------------------------------------------------------------------------

plot_title = "Nb recaptures\nas adult female from\nadult female at t+2"
sim_distribution = unlist(lapply(total_nb_recaptures_af_from_af_t2_simulated, sum))
true_value = sum(total_nb_recaptures_af_from_af_t2_true)

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.71. Total number of dead recoveries (11) ----
# ---------------------------------------------

plot_title = "Number of\ndead recoveries"
sim_distribution = total_nb_dead_recoveries_simulated
true_value = total_nb_dead_recoveries_true

sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = plot_title,
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.72. Plot ----
# -------------

# Prepare dataframe
metric_0 = names(table(sim_obs_df$metric[which(sim_obs_df$sim_distribution == 0)]))[which(table(sim_obs_df$metric[which(sim_obs_df$sim_distribution == 0)]) == 5000)]
sim_obs_df = sim_obs_df[- which(sim_obs_df$metric %in% metric_0), ]

sim_obs_df = rbind(sim_obs_df, data.frame(metric = rep(c("1\n1",
                                                         "2\n2",
                                                         "3\n3"), each = 5000),
                                          sim_distribution = rep(NA, 5000 * 3),
                                          true_value = rep(NA, 5000 * 3)))

temp = seq(1, 64, 16) # Sequence of numbers to know how many plots
                      # should appear on each page

# Factorize metric names
sim_obs_df$metric = factor(sim_obs_df$metric, levels = unique(sim_obs_df$metric))


fontSize = 7 # Smaller font size to avoid label overlap

for(i in 1:length(temp)){
  
  j = ifelse(i != length(temp), temp[i] + 15, # Index of last plot on the current page
             length(unique(sim_obs_df$metric))) 
  
  sim_obs_plot = ggplot(sim_obs_df[which(sim_obs_df$metric %in% unique(sim_obs_df$metric)[(temp[i]):(j)]), ], 
                        aes(x = sim_distribution)) +
    facet_wrap(~ metric, scales = "free", ncol = 4) +
    geom_density(alpha = 0.2, col = "#3B0F70FF", fill = "#3B0F70FF") +
    geom_vline(aes(xintercept = true_value), 
               linetype = "solid", color = "#FE9F6DFF", size = 1) +
    xlab("Parameter value") +
    ylab("Density") +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    theme_minimal() %+replace%    
    
    theme(axis.text = element_text(colour = colPlot, size = fontSize, family = font),
          axis.title.x = element_text(colour = colPlot, size = fontSize, family = font,
                                      margin = margin(t = 15, r = 0, b = 0, l = 0)), 
          axis.title.y = element_text(colour = colPlot, size = fontSize, family = font,
                                      margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = colPlot),
          axis.ticks = element_line(colour = colPlot),
          plot.background = element_rect(fill = colBG, color = colBG, size = 0),
          plot.title = element_text(colour = colPlot, size = fontSize, family = font,
                                    face = "bold", hjust = 0, 
                                    margin = margin(t = 0, r = 0, b = 15, l = 0)),
          legend.position = "none",
          strip.text = element_text(size = fontSize, margin = margin(b = 5)))
  
  png(filename = paste0("Output/Plots/Lions_Multistate_PPC_Plots_", i, ".png"),
      width = 15,
      height = 12,
      units = "cm",
      bg = "transparent",
      res = 600,
      type = "cairo")
  
  print(sim_obs_plot)
  
  dev.off()
}


## 6.2. Plotting p-values ----
# -----------------------

## 6.2.1. Total number of recaptures ----
# ----------------------------------

single_pval_df = data.frame(metric = "Total nb recaptures",
                            pvalues = total_nb_recaptures_pvalue)


## 6.2.2. Total number of recaptures at t + 1 ----
# -------------------------------------------

multi_pval_df = data.frame(metric = "Nb recaptures t+1",
                           pvalues = total_nb_recaptures_t1_pvalue)


## 6.2.3. Total number of recaptures at t + 2 ----
# -------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb recaptures t+2",
                                 pvalues = total_nb_recaptures_t2_pvalue))


## 6.2.4. Total number of recaptures as female old subadults (2) ----
# --------------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb female old subadults",
                                  pvalues = total_nb_recaptures_sa2f_pvalue))


## 6.2.5. Number of recaptures as female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb female old subadults t+1",
                                pvalues = total_nb_recaptures_sa2f_t1_pvalue))


## 6.2.6. Total number of recaptures as male old subadults (3) ----
# ------------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb male old subadults",
                                  pvalues = total_nb_recaptures_sa2m_pvalue))


## 6.2.7. Number of recaptures as male old subadults (3) at t + 1 ----
# ---------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb male old subadults t+1",
                                pvalues = total_nb_recaptures_sa2f_pvalue))


## 6.2.8. Total number of recaptures as young males (5, 6, 7, 8) ----
# --------------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb young males",
                                  pvalues = total_nb_recaptures_ym_pvalue))


## 6.2.9. Number of recaptures as young males (5, 6, 7, 8) at t + 1 ----
# -----------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young males t+1",
                                pvalues = total_nb_recaptures_ym_t1_pvalue))


## 6.2.10. Number of recaptures as young males (5, 6, 7, 8) at t + 2 ----
# ------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young males t+2",
                                  pvalues = total_nb_recaptures_ym_t2_pvalue))


## 6.2.11. Total number of recaptures as young males 1 (5) ----
# --------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb young male 1",
                                  pvalues = total_nb_recaptures_ym1_pvalue))


## 6.2.12. Number of recaptures as young males 1 (5) at t + 1 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 1 t+1",
                                  pvalues = total_nb_recaptures_ym1_t1_pvalue))


## 6.2.13. Number of recaptures as young males 1 (5) at t + 2 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 1 t+2",
                                 pvalues = total_nb_recaptures_ym1_t2_pvalue))


## 6.2.14. Number of recaptures as young males 1 (5) from male old subadult (3) at t + 1 ----
# --------------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 1 from male old subadult t+1",
                                 pvalues = total_nb_recaptures_ym1_from_sa2m_t1_pvalue))


## 6.2.15. Number of recaptures as young males 1 (5) from young subadult (1) at t + 2 ----
# -----------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 1 from young subadult t+2",
                                 pvalues = total_nb_recaptures_ym1_from_sa1_t2_pvalue))


## 6.2.16. Total number of recaptures as young males 2 (6) ----
# --------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb young male 2",
                                  pvalues = total_nb_recaptures_ym2_pvalue))


## 6.2.17. Number of recaptures as young males 2 (6) at t + 1 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 2 t+1",
                                 pvalues = total_nb_recaptures_ym2_t1_pvalue))


## 6.2.18. Number of recaptures as young males 2 (6) at t + 2 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 2 t+2",
                                 pvalues = total_nb_recaptures_ym2_t2_pvalue))


## 6.2.19. Number of recaptures as young males 2 (6) from young male 1 (5) at t + 1 ----
# ---------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 2 from young male 1 t+1",
                                 pvalues = total_nb_recaptures_ym2_from_ym1_t1_pvalue))


## 6.2.20. Number of recaptures as young males 2 (5) from male old subadult (3) at t + 2 ----
# --------------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 2 from male old subadult t+2",
                                 pvalues = total_nb_recaptures_ym2_from_sa2m_t2_pvalue))


## 6.2.21. Total number of recaptures as young males 3 (7) ----
# --------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb young male 3",
                                  pvalues = total_nb_recaptures_ym3_pvalue))


## 6.2.22. Number of recaptures as young males 3 (7) at t + 1 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 3 t+1",
                                 pvalues = total_nb_recaptures_ym3_t1_pvalue))


## 6.2.23. Number of recaptures as young males 3 (7) at t + 2 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 3 t+2",
                                 pvalues = total_nb_recaptures_ym3_t2_pvalue))


## 6.2.24. Number of recaptures as young males 3 (7) from young male 2 (6) at t + 1 ----
# ---------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 3 from young male 2 t+1",
                                 pvalues = total_nb_recaptures_ym3_from_ym2_t1_pvalue))


## 6.2.25. Number of recaptures as young males 3 (7) from young male 1 (5) at t + 2 ----
# ---------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 3 from young male 1 t+2",
                                 pvalues = total_nb_recaptures_ym3_from_ym1_t2_pvalue))


## 6.2.26. Total number of recaptures as young males 4 (8) ----
# --------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb young male 4",
                                  pvalues = total_nb_recaptures_ym4_pvalue))


## 6.2.27. Number of recaptures as young males 4 (8) at t + 2 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 4 t+2",
                                 pvalues = total_nb_recaptures_ym4_t2_pvalue))


## 6.2.28. Number of recaptures as young males 4 (8) from young male 2 (6) at t + 2 ----
# ---------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb young male 4 from young male 2 t+2",
                                 pvalues = total_nb_recaptures_ym4_from_ym2_t2_pvalue))


## 6.2.29. Total number of recaptures as nomads (9) ----
# -------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb nomadic males",
                                  pvalues = total_nb_recaptures_nm_pvalue))


## 6.2.30. Number of recaptures as nomads (9) at t + 1 ----
# ----------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic males t+1",
                                 pvalues = total_nb_recaptures_nm_t1_pvalue))


## 6.2.31. Number of recaptures as nomads (9) at t + 2 ----
# ----------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male t+2",
                                 pvalues = total_nb_recaptures_nm_t2_pvalue))


## 6.2.32. Number of recaptures as nomads (9) from male old subadult (3) at t + 2 ----
# -------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from male old subadult t+2",
                                 pvalues = total_nb_recaptures_nm_from_sa2m_t2_pvalue))


## 6.2.33. Number of recaptures as nomads (9) from young male 1 (5) at t + 1 ----
# --------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from young male 1 t+1",
                                 pvalues = total_nb_recaptures_nm_from_ym1_t1_pvalue))


## 6.2.34. Number of recaptures as nomads (9) from young male 1 (5) at t + 2 ----
# --------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from young male 1 t+2",
                                 pvalues = total_nb_recaptures_nm_from_ym1_t2_pvalue))


## 6.2.35. Number of recaptures as nomads (9) from young male 2 (6) at t + 1 ----
# --------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from young male 2 t+1",
                                 pvalues = total_nb_recaptures_nm_from_ym2_t1_pvalue))


## 6.2.36. Number of recaptures as nomads (9) from young male 2 (6) at t + 2 ----
# --------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from young male 2 t+2",
                                 pvalues = total_nb_recaptures_nm_from_ym2_t2_pvalue))


## 6.2.37. Number of recaptures as nomads (9) from nomadic male (9) at t + 1 ----
# --------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from nomadic male t+1",
                                 pvalues = total_nb_recaptures_nm_from_nm_t1_pvalue))


## 6.2.38. Number of recaptures as nomads (9) from nomadic male (9) at t + 2 ----
# --------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from nomadic male t+2",
                                 pvalues = total_nb_recaptures_nm_from_nm_t2_pvalue))


## 6.2.39. Number of recaptures as nomads (9) from resident male (10) at t + 1 ----
# ----------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from resident male t+1",
                                 pvalues = total_nb_recaptures_nm_from_rm_t1_pvalue))


## 6.2.40. Number of recaptures as nomads (9) from resident male (10) at t + 2 ----
# ----------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb nomadic male from resident male t+2",
                                 pvalues = total_nb_recaptures_nm_from_rm_t2_pvalue))


## 6.2.41. Total number of recaptures as residents (10) ----
# -----------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb resident males",
                                  pvalues = total_nb_recaptures_rm_pvalue))


## 6.2.42. Number of recaptures as residents (10) at t + 1 ----
# --------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male t+1",
                                 pvalues = total_nb_recaptures_rm_t1_pvalue))


## 6.2.43. Number of recaptures as residents (10) at t + 2 ----
# --------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male t+2",
                                 pvalues = total_nb_recaptures_rm_t2_pvalue))


## 6.2.44. Number of recaptures as residents (10) from male old subadult (3) at t + 2 ----
# -----------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from male old subadult t+2",
                                 pvalues = total_nb_recaptures_rm_from_sa2m_t2_pvalue))


## 6.2.45. Number of recaptures as residents (10) from young male 1 (5) at t + 1 ----
# ------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from young male 1 t+1",
                                 pvalues = total_nb_recaptures_rm_from_ym1_t1_pvalue))


## 6.2.46. Number of recaptures as residents (10) from young male 1 (5) at t + 2 ----
# ------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from young male 1 t+2",
                                 pvalues = total_nb_recaptures_rm_from_ym1_t2_pvalue))


## 6.2.47. Number of recaptures as residents (10) from young male 2 (6) at t + 1 ----
# ------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from young male 2 t+1",
                                 pvalues = total_nb_recaptures_rm_from_ym2_t1_pvalue))


## 6.2.48. Number of recaptures as residents (10) from young male 2 (6) at t + 2 ----
# ------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from young male 2 t+2",
                                 pvalues = total_nb_recaptures_rm_from_ym2_t2_pvalue))


## 6.2.49. Number of recaptures as residents (10) from nomadic male (9) at t + 1 ----
# ------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from nomadic male t+1",
                                 pvalues = total_nb_recaptures_rm_from_nm_t1_pvalue))


## 6.2.50. Number of recaptures as residents (10) from nomadic male (9) at t + 2 ----
# ------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from nomadic male t+2",
                                 pvalues = total_nb_recaptures_rm_from_nm_t2_pvalue))


## 6.2.51. Number of recaptures as residents (10) from resident male (10) at t + 1 ----
# --------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from resident male t+1",
                                 pvalues = total_nb_recaptures_rm_from_rm_t1_pvalue))


## 6.2.52. Number of recaptures as residents (10) from resident male (10) at t + 2 ----
# --------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb resident male from resident male t+2",
                                 pvalues = total_nb_recaptures_rm_from_rm_t2_pvalue))


## 6.2.53. Total number of recaptures as adult females (4) ----
# --------------------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb adult females",
                                  pvalues = total_nb_recaptures_af_pvalue))


## 6.2.54. Number of recaptures as adult females (4) at t + 1 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult females t+1",
                                 pvalues = total_nb_recaptures_af_t1_pvalue))


## 6.2.55. Number of recaptures as adult females (4) from female old subadults (2) at t + 1 ----
# -----------------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult female from female old subadult t+1",
                                 pvalues = total_nb_recaptures_af_from_sa2f_t1_pvalue))


## 6.2.56. Number of recaptures as adult females (4) at t + 2 ----
# -----------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult females t+2",
                                 pvalues = total_nb_recaptures_af_t2_pvalue))


## 6.2.57. Number of recaptures as adult females (4) from young subadults (1) at t + 2 ----
# ------------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult female from young subadult t+2",
                                 pvalues = total_nb_recaptures_af_from_sa1_t2_pvalue))


## 6.2.58. Number of recaptures as adult females (4) from female old subadults (2) at t + 2 ----
# -----------------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult female from female old subadult t+2",
                                 pvalues = total_nb_recaptures_af_from_sa2f_t2_pvalue))


## 6.2.59. Number of recaptures as adult females (4) from adult females (4) at t + 1 ----
# ----------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult frmale from adult female t+1",
                                 pvalues = total_nb_recaptures_af_from_af_t1_pvalue))


## 6.2.60. Number of recaptures as adult females (4) from adult females (4) at t + 2 ----
# ----------------------------------------------------------------------------------

multi_pval_df = rbind(multi_pval_df, 
                      data.frame(metric = "Nb adult female from adult female t+2",
                                 pvalues = total_nb_recaptures_af_from_af_t2_pvalue))


## 6.2.61. Total number of dead recoveries (11) ----
# ---------------------------------------------

single_pval_df = rbind(single_pval_df, 
                       data.frame(metric = "Nb dead recoveries",
                                  pvalues = total_nb_dead_recoveries_pvalue))


# Factorize the metric names
single_pval_df$metric = factor(single_pval_df$metric, levels = unique(single_pval_df$metric))
multi_pval_df$metric = factor(multi_pval_df$metric, levels = unique(multi_pval_df$metric))


## 6.3. Plot single p-values ----
# --------------------------

png(filename = "Output/Plots/Lions_Multistate_PPC_SinglePvalues.png",
    width = 8,
    height = 8,
    units = "cm",
    bg = "transparent",
    res = 600,
    type = "cairo")

ggplot(single_pval_df, aes(x = metric, y = pvalues)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = 0.5), 
             linetype = "solid", color = "#FE9F6DFF", size = 1) +
  xlab("Metric") +
  ylab("Bayesian p-value") +
  ylim(0, 1) +
  theme_minimal() %+replace%
  
  theme(axis.text = element_text(colour = colPlot, size = fontSize, family = font),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9),
        axis.title.x = element_text(colour = colPlot, size = fontSize, family = font, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = fontSize, family = font, margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, family = font, face = "bold", hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "none",
        strip.text = element_text(size = fontSize, margin = margin(b = 5)))

dev.off()


# Calculate 95% quantile for multiple p-values
mean_multi_pval_df = aggregate(pvalues ~ metric, data = multi_pval_df,
                               FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


## 6.4. Plot multi p-values ----
# -------------------------
fontSize = 6

png(filename = "Output/Plots/Lions_Multistate_PPC_MultiPvalues.png",
    width = 23,
    height = 10,
    units = "cm",
    bg = "transparent",
    res = 600,
    type = "cairo")

ggplot(multi_pval_df, aes(x = metric, y = pvalues)) +
  geom_boxplot(ymin = mean_multi_pval_df$pvalues[, 1], ymax = mean_multi_pval_df$pvalues[, 3], outlier.alpha = 0.2, size = 0.5) +
  geom_hline(aes(yintercept = 0.5), 
             linetype = "solid", color = "#FE9F6DFF", size = 1) +
  xlab("Metric") +
  ylab("Bayesian p-value") +
  ylim(0, 1) +
  theme_minimal() %+replace%    
  
  theme(axis.text = element_text(colour = colPlot, size = fontSize, family = font),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9),
        axis.title.x = element_text(colour = colPlot, size = fontSize, family = font, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = fontSize, family = font, margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, family = font, face = "bold", hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "none",
        strip.text = element_text(size = fontSize, margin = margin(b = 5)))

dev.off()
