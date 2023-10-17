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
library(viridis)


## 1.3. Loading data ----
# ------------------

# Lion demographic dataset
lions.data = read.csv("Data/01_LionsDemographicData.csv")

# MCMC samples
load("LionFullModel_Repro_Rec_Output.RData")
lions_output_fullGLMM = as.matrix(rbind(lions_output_fullGLMM[[1]],
                                        lions_output_fullGLMM[[2]],
                                        lions_output_fullGLMM[[3]],
                                        lions_output_fullGLMM[[4]]))

# Simulated datasets
load("Repro_Rec_Simulated_Data.RData")




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
# 3. Functions to calculate summary measures on datasets ----
#
###########################################################################

## 3.1. Overall probability of reproduction  ----
# -----------------------------------------

repro.prob = function(data){
  
  return(mean(data$reproduction, na.rm = T))
  
}


## 3.2. Average age of reproducing females  ----
# ----------------------------------------

age.repro.females = function(data){
  
  return(mean(data$age.at.capture[which(data$reproduction == 1)] * 12)) # Age in months
  
}


## 3.3. Mean number of cubs  ----
# -------------------------

mean.nb.cubs = function(data){
  
  return(mean(data$cubs, na.rm = T))
  
}




###########################################################################
#
# 4. Calculate summary measures on true and simulated datasets ----
#
###########################################################################

## 4.1. True data ----
# ---------------

## 4.1.1. Overall probability of reproduction  ----
# -------------------------------------------

repro_prob_true = repro.prob(females.data)


## 4.1.2. Average age of reproducing females  ----
# ------------------------------------------

age_repro_females_true = age.repro.females(females.data)


## 4.1.3. Mean number of cubs  ----
# ---------------------------

mean_nb_cubs_true = mean.nb.cubs(females.data)


## 4.2. Simulated data ----
# --------------------

## 4.2.1. Overall probability of reproduction  ----
# -------------------------------------------

repro_prob_simulated = c()

for(nsim in 1:length(list_simulated_data_seasonal_coeffs)){
  
  repro_prob_simulated = c(repro_prob_simulated, 
                           repro.prob(list_simulated_data_seasonal_coeffs[[nsim]]))
  
} 


## 4.2.2. Average age of reproducing females  ----
# ------------------------------------------

age_repro_females_simulated = c()

for(nsim in 1:length(list_simulated_data_seasonal_coeffs)){
  
  age_repro_females_simulated = c(age_repro_females_simulated, 
                                  age.repro.females(list_simulated_data_seasonal_coeffs[[nsim]]))
  
} 


## 4.2.3. Mean number of cubs  ----
# ---------------------------

mean_nb_cubs_simulated = c()

for(nsim in 1:length(list_simulated_data_seasonal_coeffs)){
  
  mean_nb_cubs_simulated = c(mean_nb_cubs_simulated, 
                             mean.nb.cubs(list_simulated_data_seasonal_coeffs[[nsim]]))
  
}




###########################################################################
#
# 5. Calculating Bayesian p-values ----
#
###########################################################################

## 5.1. Bayesian p-values ----
# -----------------------

## 5.1.1. Overall probability of reproduction  ----
# -------------------------------------------

repro_prob_pvalue = length(which(repro_prob_simulated > repro_prob_true)) / 
                    length(list_simulated_data_seasonal_coeffs)


## 5.1.2. Average age of reproducing females  ----
# ------------------------------------------

age_repro_females_pvalue = length(which(age_repro_females_simulated > age_repro_females_true)) /
                           length(list_simulated_data_seasonal_coeffs)


## 5.1.3. Mean number of cubs  ----
# ---------------------------

mean_nb_cubs_pvalue = length(which(mean_nb_cubs_simulated > mean_nb_cubs_true)) / 
                      length(list_simulated_data_seasonal_coeffs)




###########################################################################
#
# 6. Plotting simulated and true metric values and Bayesian p-values ----
#
###########################################################################

## 6.1. Plotting simulated and true value ----
# ---------------------------------------

# Plot settings :

# Define background color of plot 
colBG = "transparent"

# Define color, font, font size of plot
colPlot = "black"
font = "Helvetica"
fontSize = 10


# We first merge all the simulated and observed values in the same
# dataframe before plotting everything together.


## 6.1.1. Overall probability of reproduction  ----
# -------------------------------------------

plot_title = "True vs. simulated value: Proportion of females reproducing"
sim_distribution = repro_prob_simulated
true_value = repro_prob_true

prior_post = data.frame(metric = "Proportion of\nfemales reproducing",
                        sim_distribution = sim_distribution,
                        true_value = true_value)


## 6.1.2. Average age of reproducing females  ----
# ------------------------------------------

plot_title = "True vs. simulated value: Mean age of reproducing females (in months)"
sim_distribution = age_repro_females_simulated
true_value = age_repro_females_true

prior_post = rbind(prior_post, 
                   data.frame(metric = "Mean age of\nreproducing females\n(in months)",
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.3. Mean number of cubs  ----
# ---------------------------

plot_title = "True vs. simulated value: Mean number of cubs per female"
sim_distribution = mean_nb_cubs_simulated
true_value = mean_nb_cubs_true

prior_post = rbind(prior_post,
                   data.frame(metric = "Mean number of\ncubs per female",
                              sim_distribution = sim_distribution,
                              true_value = true_value))


## 6.1.4. Plot ----
# ------------

# Prepare dataframe
sim_obs_df = rbind(sim_obs_df,
                   data.frame(metric = rep(c("Nb recaptures\nas resident male from\nnomadic male at t+1", 
                                             "Nb recaptures\nas resident male from\nnomadic male at t+2", 
                                             "Nb recaptures\nas resident male from\nresident male at t+1",
                                             "Nb recaptures\nas resident male from\nresident male at t+2",
                                             "Nb recaptures\nas adult female",                            
                                             "Nb recaptures as\nadult female at t+1",                     
                                             "Nb recaptures as\nadult female from\n old subadult\nat t+1",
                                             "Nb recaptures as\nadult female at t+2",                     
                                             "Nb recaptures\nas adult female from\nyoung subadult at t+2",
                                             "Nb recaptures as\nadult female from\n old subadult\nat t+2",
                                             "Nb recaptures\nas adult female from\nadult female at t+1",  
                                             "Nb recaptures\nas adult female from\nadult female at t+2",  
                                             "Number of\ndead recoveries"), each = 5000),
                              sim_distribution = rep(NA, 5000 * 13),
                              true_value = rep(NA, 5000 * 13)))

# Factorize metric names
sim_obs_df$metric = factor(sim_obs_df$metric, unique(sim_obs_df$metric))


png(filename = "PPC_Plots.png", 
    width = 15, 
    height = 12, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(sim_obs_df, aes(x = sim_distribution)) +
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
        axis.title.x = element_text(colour = colPlot, size = fontSize, 
                                    family = font, 
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = fontSize, 
                                    family = font, 
                                    margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, 
                                  family = font, face = "bold", hjust = 0,
                                  margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "none",
        strip.text = element_text(size = fontSize, margin = margin(b = 5)))

dev.off()


## 6.2. Plotting p-values ----
# -----------------------

## 6.2.1. Overall probability of reproduction ----
# -------------------------------------------

pval_df = data.frame(metric = "Proportion\n of females\nreproducing",
                     pvalues = repro_prob_pvalue)


## 6.2.2. Average age of reproducing females ----
# ------------------------------------------

pval_df = rbind(pval_df, 
                data.frame(metric = "Mean age of\nreproducing\nfemales",
                           pvalues = age_repro_females_pvalue))


## 6.2.3. Mean number of cubs ----
# ---------------------------

pval_df = rbind(pval_df, 
                data.frame(metric = "Mean number of\ncubs per female",
                           pvalues = mean_nb_cubs_pvalue))


## 6.3. Plot p-values ----
# -------------------

png(filename = "ReproRec_PPC_Pvalues.png",
    width = 8,
    height = 8,
    units = "cm",
    bg = "transparent",
    res = 600,
    type = "cairo")

ggplot(pval_df, aes(x = metric, y = pvalues)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = 0.5), linetype = "solid", 
             color = "#FE9F6DFF", size = 1) +
  xlab("Metric") +
  ylab("Bayesian p-value") +
  ylim(0, 1) +
  theme_minimal() %+replace%    
  theme(axis.text = element_text(colour = colPlot, size = fontSize, family = font),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9),
        axis.title.x = element_text(colour = colPlot, size = fontSize, 
                                    family = font, 
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(colour = colPlot, size = fontSize,
                                    family = font, 
                                    margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
        plot.title = element_text(colour = colPlot, size = fontSize, 
                                  family = font, face = "bold", hjust = 0, 
                                  margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "none",
        strip.text = element_text(size = fontSize, margin = margin(b = 5)))

dev.off()




###########################################################################
#
# 7. Plotting observed and predicted values ----
#
###########################################################################

## 7.1. Function to predict reproduction probability and recruitment ----
# ------------------------------------------------------------------

predict_repro_rec = function(output, # Vital rate to predict (repro or rec)
                             
                             # Intercepts of vital-rate models
                             mu.repro = lions_output_fullGLMM[, grep("mu.repro", # Reproduction probability
                                                                     colnames(lions_output_fullGLMM))],
                             mu.rec = lions_output_fullGLMM[, grep("mu.rec",     # Recruitment
                                                                   colnames(lions_output_fullGLMM))],
                             
                             # Betas of vital-rate models
                             # Reproduction probability
                             repro.beta.habitat.woodland = lions_output_fullGLMM[, grep("repro.beta.habitat.woodland",                   # Habitat
                                                                                        colnames(lions_output_fullGLMM))],
                             repro.beta.age = lions_output_fullGLMM[, grep("repro.beta.age",                                             # Age
                                                                           colnames(lions_output_fullGLMM))],
                             repro.beta.quad.age = lions_output_fullGLMM[, grep("repro.beta.quad.age",                                   # Quadratic age
                                                                                colnames(lions_output_fullGLMM))],
                             repro.beta.nb.af.pride = lions_output_fullGLMM[, grep("repro.beta.nb.af.pride",                             # Number of adult females in the pride
                                                                                   colnames(lions_output_fullGLMM))][, c(1, 2)],
                             repro.beta.quad.nb.af.pride = lions_output_fullGLMM[, grep("repro.beta.quad.nb.af.pride",                   # Quadratic number of adult females in the pride
                                                                                        colnames(lions_output_fullGLMM))],
                             repro.beta.nb.nm.coal.hr = lions_output_fullGLMM[, grep("repro.beta.nb.nm.coal.hr",                         # Number of nomadic coalitions in the home range
                                                                                     colnames(lions_output_fullGLMM))], 
                             repro.beta.nb.af.pride.nb.nm.coal.hr = lions_output_fullGLMM[, grep("repro.beta.nb.af.pride.nb.nm.coal.hr", # Interaction between the number of
                                                                                                                                         # adult females in the pride and the
                                                                                                                                         # number of nomadic coalitions in 
                                                                                                                                         # the home range 
                                                                                                 colnames(lions_output_fullGLMM))],
                             # Recruitment
                             rec.beta.habitat.woodland = lions_output_fullGLMM[, grep("rec.beta.habitat.woodland",                       # Habitat
                                                                                      colnames(lions_output_fullGLMM))],
                             rec.beta.nb.af.pride = lions_output_fullGLMM[, grep("rec.beta.nb.af.pride",                                 # Number of adult females in the pride
                                                                                 colnames(lions_output_fullGLMM))][, c(1, 2)],
                             rec.beta.nb.nm.coal.hr = lions_output_fullGLMM[, grep("rec.beta.nb.nm.coal.hr",                             # Number of nomadic coalitions in the home range
                                                                                   colnames(lions_output_fullGLMM))],
                             rec.beta.nb.af.pride.nb.nm.coal.hr = lions_output_fullGLMM[, grep("rec.beta.nb.af.pride.nb.nm.coal.hr",     # Interaction between the number of
                                                                                                                                         # adult females in the pride and the
                                                                                                                                         # number of nomadic coalitions in 
                                                                                                                                         # the home range 
                                                                                               colnames(lions_output_fullGLMM))],
                             
                             # Random effect epsilons
                             epsilon.repro = lions_output_fullGLMM[, grep("epsilon.repro", # Reproduction probability
                                                                          colnames(lions_output_fullGLMM))],
                             epsilon.rec = lions_output_fullGLMM[, grep("epsilon.rec",     # Recruitment
                                                                        colnames(lions_output_fullGLMM))],
                             
                             # Covariates
                             season,        # Season at a given timestep
                             year,          # Year at a given timestep
                             nb.nm.coal.hr, # Number of nomadic coalitions in the home range
                             nb.af.pride,   # Number of adult females in the pride
                             habitat,       # Habitat
                             habitat.prob = habitat.prob, # Habitat probability for missing habitat values
                             age){          # Age
  
  # Sample habitat missing value
  habitat = ifelse(!is.na(habitat), habitat, 
                   rbinom(1, size = 1, prob = habitat.prob) + 1)
  
  # Initialize vector of vital rate
  repro = rec = c()
  
  
  # Calculate reproduction probability prediction
  if(output == "repro"){
    
    for(i in 1:nrow(mu.repro)){
      
      repro = c(repro, as.numeric(plogis(mu.repro[i, season] + 
                                         matrix(repro.beta.habitat.woodland[i, ], 2, 2)[season, habitat] +
                                         repro.beta.age[i, season] * age +
                                         repro.beta.quad.age[i, season] * age^2 +
                                         repro.beta.nb.af.pride[i, season] * nb.af.pride +
                                         repro.beta.quad.nb.af.pride[i, season] * nb.af.pride^2 +
                                         repro.beta.nb.nm.coal.hr[i, season] * nb.nm.coal.hr +
                                         repro.beta.nb.af.pride.nb.nm.coal.hr[i, season] * nb.af.pride * nb.nm.coal.hr +
                                         epsilon.repro[i, grep(paste("epsilon.repro\\[", season, ", ", year, "\\]", sep = ""), colnames(epsilon.repro))])))
    }
    
    return(repro)
  }
  
  
  # Calculate recruitment prediction
  else if(output == "rec"){
    
    for(i in 1:nrow(mu.repro)){
      
      rec = c(rec, as.numeric(exp(mu.rec[i, season] + 
                                  matrix(rec.beta.habitat.woodland[i, ], 2, 2)[season, habitat] +
                                  rec.beta.nb.af.pride[i, season] * nb.af.pride +
                                  rec.beta.nb.nm.coal.hr[i, season] * nb.nm.coal.hr +
                                  rec.beta.nb.af.pride.nb.nm.coal.hr[i, season] * nb.af.pride * nb.nm.coal.hr +
                                  epsilon.rec[i, grep(paste("epsilon.rec\\[", season, ", ", year, "\\]", sep = ""), colnames(epsilon.rec))])))
      
    }
    
    return(rec)
  }
}


## 7.2. Reproduction probability ----
# ------------------------------

# Observed reproduction probability
repro_obs = aggregate(reproduction ~ year + season, 
                      data = females.data, mean)


# Initialize reproduction probability prediction dataframe
repro_pred = repro_obs[, 1:2]
repro_pred$lwr = NA
repro_pred$mean = NA
repro_pred$upr = NA

# Fill in reproduction probability prediction mean and credible intervals
for(i in 1:nrow(repro_pred)){

  print(i)
  
  # Get season and year
  seas = unique(females.data$season.nb[which(females.data$season == repro_pred$season[i])])
  yr = unique(females.data$year.nb[which(females.data$year == repro_pred$year[i])])

  # Subset data to current season and year
  females_subset = females.data[which(females.data$season.nb == seas & females.data$year.nb == yr), ]

  # Initialize reproduction probability vector
  repro = c()
  
  # Predict reproduction for each row of the subset
  for(j in 1:nrow(females_subset)){

    repro = c(repro, predict_repro_rec(output = "repro",
                                       season = seas,
                                       year = yr,
                                       nb.nm.coal.hr = females_subset$nb.nm.coal.hr.scaled[j],
                                       nb.af.pride = females_subset$nb.af.pride.scaled[j],
                                       habitat = females_subset$habitat.code[j],
                                       habitat.prob = habitat.prob,
                                       age = females_subset$age.at.capture.scaled[j]))

  }

  # Get mean and 95% quantile of prediction
  repro_pred$mean[i] = mean(repro, na.rm = T)
  repro_pred$lwr[i] = quantile(repro, probs = 0.025, na.rm = T)
  repro_pred$upr[i] = quantile(repro, probs = 0.975, na.rm = T)
}

repro_pred

# Save predictions
save(repro_pred, file = "Repro_Pred.RData")


## 7.3. Recruitment ----
# -----------------

# Observed recruitment
rec_obs = aggregate(cubs ~ year + season, 
                    data = females.data, mean)


# Initialize recruitment prediction dataframe
rec_pred = rec_obs[, 1:2]
rec_pred$lwr = NA
rec_pred$mean = NA
rec_pred$upr = NA

# Fill in recruitment prediction mean and credible intervals
for(i in 1:nrow(rec_pred)){
  
  print(i)
  
  # Get season and year
  seas = unique(females.data$season.nb[which(females.data$season == rec_pred$season[i])])
  yr = unique(females.data$year.nb[which(females.data$year == rec_pred$year[i])])
  
  # Subset data to current season and year
  females_subset = females.data[which(females.data$season.nb == seas & females.data$year.nb == yr), ]
  
  # Initialize recruitment vector
  rec = c()
  
  # Predict recruitment for each row of the subset
  for(j in 1:nrow(females_subset)){
    
    rec = c(rec, predict_repro_rec(output = "rec",
                                   season = seas,
                                   year = yr,
                                   nb.nm.coal.hr = females_subset$nb.nm.coal.hr.scaled[j],
                                   nb.af.pride = females_subset$nb.af.pride.scaled[j],
                                   habitat = females_subset$habitat.code[j],
                                   habitat.prob = habitat.prob,
                                   age = females_subset$age.at.capture.scaled[j]))
    
  }
  
  # Get mean and 95% quantile of prediction
  rec_pred$mean[i] = mean(rec, na.rm = T)
  rec_pred$lwr[i] = quantile(rec, probs = 0.025, na.rm = T)
  rec_pred$upr[i] = quantile(rec, probs = 0.975, na.rm = T)
}

rec_pred

# Save predictions
save(rec_pred, file = "Rec_Pred.RData")


## 7.4. Plot observed vs predicted values ----
# ---------------------------------------

# Add index to data for plotting
repro_pred$timestep = paste(repro_pred$year, repro_pred$season, sep = "_")
repro_pred$index[order(repro_pred$timestep)] = seq(1, nrow(repro_pred))

rec_pred$timestep = paste(rec_pred$year, rec_pred$season, sep = "_")
rec_pred$index[order(rec_pred$timestep)] = seq(1, nrow(rec_pred))


# Colourblind palette
cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Plot observed vs predicted reproduction probabilities
png(file = "ObsVSPred_Repro.png",
    type = "cairo",
    units = "cm",
    width = 15,
    height = 10,
    res = 600,
    bg = "transparent")

ggplot(repro_pred, aes(x = year, y = mean, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2, col = NA) +
  geom_point(aes(x = repro_obs$year, y = repro_obs$reproduction, 
                 colour = repro_obs$season), size = 2) +
  scale_color_manual(values = c(cbbPalette[2], cbbPalette[4]),
                     label = c("Dry", "Wet"),
                     name = "Season") +
  scale_fill_manual(values = c(cbbPalette[2], cbbPalette[4]),
                    label = c("Dry", "Wet"),
                    name = "Season") +
  xlab("Year") +
  ylab("Proportion of females reproducing") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, colour = "black", family = font),
        axis.title = element_text(size = 10, family = font),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        legend.title = element_text(size = 10, family = font),
        legend.text = element_text(size = 10, family = font),
        axis.title.x = element_text(vjust = -0.4),
        axis.title.y = element_text(vjust = 0.9),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))

dev.off()


# Plot observed vs predicted recruitment
png(file = "ObsVSPred_Rec.png",
    type = "cairo",
    units = "cm",
    width = 15,
    height = 10,
    res = 600,
    bg = "transparent")

ggplot(rec_pred, aes(x = year, y = mean, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season),
              alpha = 0.2, col = NA) +
  geom_point(aes(x = rec_obs$year, y = rec_obs$cubs, colour = rec_obs$season), 
             size = 2) +
  scale_color_manual(values = c(cbbPalette[2], cbbPalette[4]),
                     label = c("Dry", "Wet"),
                     name = "Season") +
  scale_fill_manual(values = c(cbbPalette[2], cbbPalette[4]),
                    label = c("Dry", "Wet"),
                    name = "Season") +
  xlab("Year") +
  ylab("Recruitment to 1 year old") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, colour = "black", family = font),
        axis.title = element_text(size = 10, family = font),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        legend.title = element_text(size = 10, family = font),
        legend.text = element_text(size = 10, family = font),
        axis.title.x = element_text(vjust = -0.4),
        axis.title.y = element_text(vjust = 0.9),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))

dev.off()
