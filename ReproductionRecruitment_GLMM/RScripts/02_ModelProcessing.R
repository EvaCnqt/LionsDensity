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

rm(list = ls(all.names = T))

## 1.1. Loading libraries ----
# -----------------------

library(ggplot2)
library(MCMCvis)
library(coda)
library(bayestestR)
library(Cairo)


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

# Format age and habitat
lions.data$age.at.capture = lions.data$age.at.capture / 12 # Age in years
lions.data$habitat.code = lions.data$habitat.code + 1 # Habitat as c(1, 2) instead of c(0, 1)


# Subset female data only and remove nomadic females
females.data = lions.data[which(lions.data$stage == "AF"), ]
females.data = females.data[- which(females.data$pride == "NO"), ]


# New season number column
females.data$season.nb = NA
females.data$season.nb[which(females.data$season == "wet")] = 1
females.data$season.nb[which(females.data$season == "dry")] = 2


# New year number column
year.number = data.frame(cbind(unique(females.data$year), # Assigning a number to each year of the dataset
                               seq(1:length(unique(females.data$year))))) 
colnames(year.number) = c("year", "year.nb")


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
# 3. Checking for convergence by plotting all three chains ----
# for each parameter
#
###########################################################################

# List of parameters
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


# Remove epsilons and sigmas for the traceplots
params.to.plot = parameters[- grep("sigma", parameters)]
params.to.plot = params.to.plot[- grep("epsilon", params.to.plot)]


# Traceplots of all chains
MCMCtrace(lions_results_fullGLMM, 
          priors = runif(45000*4, -10, 10), # Priors to get the prior-posterior overlaps 
          params = params.to.plot[grep("repro", params.to.plot)], # Parameters to plot
          iter = 45000, # Plot all iterations 
          filename = "Repro_MCMCtrace.pdf")

MCMCtrace(lions_results_fullGLMM, 
          priors = runif(45000*4, -5, 1), # Priors to get the prior-posterior overlaps 
          params = params.to.plot[grep("rec", params.to.plot)], # Parameters to plot
          iter = 45000, # Plot all iterations 
          filename = "Rec_MCMCtrace.pdf")




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
          legend.key.width = unit(35, "pt")) 
}




###########################################################################
#
# 5. Prior-posterior overlap ----
#
###########################################################################

params.prior.post.overlap = parameters[- grep("sigma", parameters)]
params.prior.post.overlap = params.prior.post.overlap[- grep("epsilon", params.prior.post.overlap)]

# Empty dataframe to store the prior and posterior data
prior_post_param_df = data.frame(parameter = NA,
                                 season = NA,
                                 prior_post = NA,
                                 distribution = NA)[0, ]


# Parameters name to plot
param_full_name = data.frame(param_short = params.prior.post.overlap,
                             param_full = c("Mean reproduction\nprobability",
                                            "Mean recruitment",
                                            "Reproduction\nHabitat (woodland)",
                                            "Reproduction\nAge",
                                            "Reproduction\nAge^2",
                                            "Reproduction\nNb females pride",
                                            "Reproduction\nNb females pride^2",
                                            "Reproduction\nNb nomadic coalitions HR",
                                            "Reproduction\nNb females\npride:Nb nomad coal. HR",
                                            "Recruitment\nHabitat (woodland)",
                                            "Recruitment\nNb females pride",
                                            "Recruitment\nNb nomadic coalitions HR",
                                            "Recruitment\nNb females\npride:Nb nomad coal. HR"))


# Fill in dataframe with prior and posterior distributions
for(param in params.prior.post.overlap){
  
  # Reproduction probability priors
  if(strsplit(param, split = ".", fixed = T)[[1]][1] == "repro" | 
     strsplit(param, split = ".", fixed = T)[[1]][2] == "repro"){
    
    dist.lwr = -10
    dist.upr = 10
  }
  
  # Recruitment priors
  else if(strsplit(param, split = ".", fixed = T)[[1]][1] == "rec" | 
          strsplit(param, split = ".", fixed = T)[[1]][2] == "rec"){
    
    dist.lwr = -5
    dist.upr = 1
  }
    
  # Full parameter name
  full_param = param_full_name$param_full[which(param_full_name$param_short == param)]
  
  # Fill in prior-posterior dataframe
  prior_post_param = data.frame(parameter = full_param,
                                season = rep(c("wet", "dry"), 
                                             each = 2 * nrow(lions_output_fullGLMM)),
                                prior_post = rep(rep(c("prior", "posterior"), 
                                                     each = nrow(lions_output_fullGLMM)), 2),
                                distribution = c(runif(nrow(lions_output_fullGLMM), dist.lwr, dist.upr), 
                                                 lions_output_fullGLMM[, grep(param, colnames(lions_output_fullGLMM))][, 1], 
                                                 runif(nrow(lions_output_fullGLMM), dist.lwr, dist.upr), 
                                                 lions_output_fullGLMM[, grep(param, colnames(lions_output_fullGLMM))][, 2]))
  
  # Merge data to the full dataframe
  prior_post_param_df = rbind(prior_post_param_df, prior_post_param)
}


# Format prior-posterior dataframe
prior_post_param_df$season = paste0("\n(", prior_post_param_df$season, " season)")
prior_post_param_df$param_seas = paste(prior_post_param_df$parameter, 
                                       prior_post_param_df$season, sep = " ")
prior_post_param_df$param_seas = factor(prior_post_param_df$param_seas, 
                                        levels = unique(prior_post_param_df$param_seas))


# Prior-posterior overlap plots
# (the percentages have been aded manually using the MCMCtrace plot information)

temp = seq(1, 13, 6) # Limits of numbers of plots per page

for(i in 1:(length(temp) - 1)){
  
  j = ifelse(i == 1, temp[i] + 5, length(unique(prior_post_param_df$parameter))) # Index of last plot on the current page
  
  # Prior-posterior overlap plot
  prior_post_plot = ggplot(prior_post_param_df[which(prior_post_param_df$parameter %in% 
                                                     unique(prior_post_param_df$parameter)[(temp[i]):(j)]), ], 
                           aes(x = distribution, fill = prior_post, colour = prior_post)) +
    facet_wrap(~ param_seas, ncol = 3, scale = "free") +
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
          plot.title = element_text(colour = colPlot, size = 8, family = font, hjust = 0, margin = margin(t = 0, r = 0, b = 15, l = 0)),
          legend.position = "right", 
          legend.justification = "right", 
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.text = element_text(colour = colPlot, size = 7, family = font),
          legend.title = element_text(colour = colPlot, size = 8, family = font),
          legend.key.size = unit(15, "pt"),
          strip.text = element_text(size = 8, margin = margin(b = 10)))
  
  if(i == 1){
    
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
  
  else{
    
    png(filename = paste0("PriorPosteriorOverlap_", i, ".png"), 
        width = 16, 
        height = 16, 
        units = "cm", 
        bg = "transparent", 
        res = 600, 
        type = "cairo")
    
    print(prior_post_plot)
    
    dev.off()
  }
}




###########################################################################
#
# 6. Effect size and prediction plots ----
#
###########################################################################

## 6.1. Reproduction probability ----
# ------------------------------

# Subset model output data and remove epsilons and sigmas
lions_output_repro = lions_output_fullGLMM[, grep("repro", colnames(lions_output_fullGLMM))]
lions_output_repro.epsilons = lions_output_repro[, grep("epsilon", colnames(lions_output_repro))]
lions_output_repro = lions_output_repro[, - grep("epsilon", colnames(lions_output_repro))]
lions_output_repro = lions_output_repro[, - grep("sigma", colnames(lions_output_repro))]


# Mean, median, and 90% credible interval of the estimates
lions_output_repro_estimates = data.frame(parameter = colnames(lions_output_repro), 
                                          mean = apply(lions_output_repro, 2, mean),
                                          median = apply(lions_output_repro, 2, median),
                                          lowerCI = apply(lions_output_repro, 2, 
                                                          FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                          upperCI = apply(lions_output_repro, 2, 
                                                          FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_repro_estimates$parameter_clean = c("Mean reproduction probability\n(wet season in the grassland)",
                                                 "Mean reproduction probability\n(dry season in the grassland)",
                                                 "Age\n(wet season)",
                                                 "Age\n(dry season)",
                                                 "Habitat\n(wet season in the grassland)",
                                                 "Habitat\n(dry season in the grassland)",
                                                 "Habitat\n(wet season in the woodland)",
                                                 "Habitat\n(dry seasons in the woodland)",
                                                 "Number of females\nin the pride (wet season)",
                                                 "Number of females\nin the pride (dry season)",
                                                 "Nb females pride:Nb nomadic\ncoal HR (wet season)",
                                                 "Nb females pride:Nb nomadic\ncoal HR (dry season)",
                                                 "Number of nomadic coalitions\nin the home range (wet season)", 
                                                 "Number of nomadic coalitions\nin the home range (dry season)",
                                                 "Age^2\n(wet season)",
                                                 "Age^2\n(dry season)",
                                                 "Nb females pride^2\n(wet season)",
                                                 "Nb females pride^2\n(dry season)")


# Add season and parameter category
lions_output_repro_estimates$season = rep(c("Wet", "Dry"), 9)
lions_output_repro_estimates$parameter_plot = rep(c("Mean reproduction probability", 
                                                    "Age", 
                                                    "Habitat (grassland)",
                                                    "Habitat (woodland)", 
                                                    "Number of females\nin the pride", 
                                                    "Nb females pride:Nb nomadic coals HR",
                                                    "Number of nomadic\ncoalitions in the home range", 
                                                    "Age^2", 
                                                    "Nb females pride^2"), each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_fullGLMM_seasonal_coeffs) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_repro_estimates$parameter_plot))),
                               variable = rep(c("Mean reproduction probability", 
                                                "Age", 
                                                "Habitat (grassland)",
                                                "Habitat (woodland)", 
                                                "Number of females\nin the pride",
                                                "Nb females pride:Nb nomadic coals HR",
                                                "Number of nomadic\ncoalitions in the home range", 
                                                "Age^2", 
                                                "Nb females pride^2"), 
                                              each = 2 * n),
                               posterior = c(lions_output_repro[1:nrow(lions_output_repro), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean reproduction probability",
                                               "Number of females\nin the pride", 
                                               "Nb females pride^2",
                                               "Number of nomadic\ncoalitions in the home range",
                                               "Nb females pride:Nb nomadic coals HR",
                                               "Habitat (woodland)", 
                                               "Age",
                                               "Age^2")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season, 
                                  levels = unique(df.plot.posterior$season))


## 6.1.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_repro_estimates$parameter_plot,
                     season = lions_output_repro_estimates$season,
                     mean = lions_output_repro_estimates$mean,
                     median = lions_output_repro_estimates$median,
                     BCI90_upper = lions_output_repro_estimates$upperCI,
                     BCI90_lower = lions_output_repro_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable,
                          levels = c("Mean reproduction probability", 
                                     "Number of females\nin the pride",
                                     "Nb females pride^2",
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Nb females pride:Nb nomadic coals HR", 
                                     "Habitat (woodland)",
                                     "Age",
                                     "Age^2")[seq(length(df.plot$variable), 1, -1)])
df.plot$season = factor(df.plot$season, 
                        levels = unique(df.plot$season))

# Set labels for plots
repro_y_labels = c("Mean reproduction probability", 
                   "Number of females\nin the pride", 
                   expression("Nb females pride"^2), 
                   "Number of nomadic\ncoalitions in the home range", 
                   "Nb females pride:Nb nomadic coals HR", 
                   "Habitat (woodland)", 
                   "Age", expression("Age"^2))
repro_y_labels = repro_y_labels[seq(length(repro_y_labels), 1, -1)]


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "Reproduction_Probability.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

repro_plot = ggplot() +
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
  scale_y_discrete(labels = repro_y_labels) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-3, 4, 1), limits = c(-4, 2)) +
  theme_general() +
  ggtitle("Reproduction probability") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

repro_plot

dev.off()


## 6.1.2. Plotting the effect sizes ----
# ---------------------------------

# Function to calculate the predictions
repro.estimate = function(season = "wet",
                          habitat = "grassland",
                          age = 0,
                          nb.af.pride = 0,
                          nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.repro = lions_output_repro[, grep("mu.repro", 
                                           colnames(lions_output_repro))][, 1]
    beta.age = lions_output_repro[, grep("repro.beta.age", 
                                         colnames(lions_output_repro))][, 1]
    beta.quad.age = lions_output_repro[, grep("repro.beta.quad.age", 
                                              colnames(lions_output_repro))][, 1]
    beta.nb.af.pride = lions_output_repro[, grep("repro.beta.nb.af.pride", 
                                                 colnames(lions_output_repro))][, 1]
    beta.quad.nb.af.pride = lions_output_repro[, grep("repro.beta.quad.nb.af.pride",
                                                      colnames(lions_output_repro))][, 1]
    beta.nb.nm.coal.hr = lions_output_repro[, grep("repro.beta.nb.nm.coal.hr", 
                                                   colnames(lions_output_repro))][, 1]
    beta.nb.af.pride.nb.nm.coal.hr = lions_output_repro[, grep("repro.beta.nb.af.pride.nb.nm.coal.hr", 
                                                               colnames(lions_output_repro))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_repro.epsilons[, seq(1, 59, 2)], 1, median) 
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_repro))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_repro[, grep("repro.beta.habitat.woodland", 
                                               colnames(lions_output_repro))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.repro = lions_output_repro[, grep("mu.repro", 
                                           colnames(lions_output_repro))][, 2]
    beta.age = lions_output_repro[, grep("repro.beta.age", 
                                         colnames(lions_output_repro))][, 2]
    beta.quad.age = lions_output_repro[, grep("repro.beta.quad.age", 
                                              colnames(lions_output_repro))][, 2]
    beta.nb.af.pride = lions_output_repro[, grep("repro.beta.nb.af.pride", 
                                                 colnames(lions_output_repro))][, 2]
    beta.quad.nb.af.pride = lions_output_repro[, grep("repro.beta.quad.nb.af.pride", 
                                                      colnames(lions_output_repro))][, 2]
    beta.nb.nm.coal.hr = lions_output_repro[, grep("repro.beta.nb.nm.coal.hr",
                                                   colnames(lions_output_repro))][, 2]
    beta.nb.af.pride.nb.nm.coal.hr = lions_output_repro[, grep("repro.beta.nb.af.pride.nb.nm.coal.hr", 
                                                               colnames(lions_output_repro))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_repro.epsilons[, seq(2, 60, 2)], 1, median) 
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_repro))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_repro[, grep("repro.beta.habitat.woodland", colnames(lions_output_repro))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred =  mean.repro + beta.habitat + beta.age * age + beta.quad.age * age^2 +
          beta.nb.af.pride * nb.af.pride + beta.quad.nb.af.pride * nb.af.pride^2 + 
          beta.nb.nm.coal.hr * nb.nm.coal.hr + 
          beta.nb.af.pride.nb.nm.coal.hr * nb.af.pride * nb.nm.coal.hr + 
          median.epsilons
  
  return(pred)
  
}


# Average reproduction probability per season and habitat

# Empty prediction dataframe
repro.season.estimates = expand.grid(season = c("wet", "dry"), 
                                     habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
repro.season.estimates$pred = exp(apply(apply(repro.season.estimates, 
                                              1, 
                                              FUN = function(x) repro.estimate(season = x[1], habitat = x[2])), 
                                        2, median))
repro.season.estimates$lwr = exp(apply(apply(repro.season.estimates, 
                                             1, 
                                             FUN = function(x) repro.estimate(season = x[1], habitat = x[2])), 
                                       2, FUN = function(x) quantile(x, probs = 0.05)))
repro.season.estimates$upr = exp(apply(apply(repro.season.estimates,
                                             1,
                                             FUN = function(x) repro.estimate(season = x[1], habitat = x[2])), 
                                       2, FUN = function(x) quantile(x, probs = 0.95)))

png(filename = "Predictions_ReproProb_SeasonHabitat.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(repro.season.estimates, aes(x = season, y = pred, colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = position_dodge(0.8), width = 0.5) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  ylim(0, 1) +
  xlab("Season") +
  ylab("Reproduction probability") +
  labs(title = "Reproduction probability") +
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


# Average reproduction probability with age

# Empty prediction dataframe
repro.age.estimates = expand.grid(season = c("wet", "dry"),
                                  age = (seq(2, 15) - mean(females.data$age.at.capture, na.rm = T)) / (2 * sd(females.data$age.at.capture, na.rm = T)))

# Fill in prediction and credible intervals
repro.age.estimates$pred = plogis(apply(apply(repro.age.estimates, 
                                              1,
                                              FUN = function(x) repro.estimate(season = x[1],
                                                                               age = rep(as.numeric(x[2]), 
                                                                                         nrow(lions_output_repro)))),
                                        2, median))
repro.age.estimates$lwr = plogis(apply(apply(repro.age.estimates, 
                                             1,
                                             FUN = function(x) repro.estimate(season = x[1], 
                                                                              age = rep(as.numeric(x[2]),
                                                                                        nrow(lions_output_repro)))), 
                                       2, FUN = function(x) quantile(x, probs = 0.05)))
repro.age.estimates$upr = plogis(apply(apply(repro.age.estimates, 
                                             1,
                                             FUN = function(x) repro.estimate(season = x[1], 
                                                                              age = rep(as.numeric(x[2]), 
                                                                                        nrow(lions_output_repro)))), 
                                       2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_ReproProb_Age.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(repro.age.estimates, aes(x = age, y = pred, colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = seq(2, 15, 2),
                     labels = seq(2, 15, 2)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Age") +
  ylab("Reproduction probability") +
  labs(title = "Reproduction probability") +
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


# Average reproduction probability with number of females in the pride

# Empty prediction dataframe
repro.nb.af.pride.estimates = expand.grid(season = c("wet", "dry"),
                                          nb.af.pride = (seq(2, 16) - mean(females.data$nb_af_pride, na.rm = T)) / (2 * sd(females.data$nb_af_pride, na.rm = T)))

# Fill in prediction and credible intervals
repro.nb.af.pride.estimates$pred = plogis(apply(apply(repro.nb.af.pride.estimates, 
                                                      1,
                                                      FUN = function(x) repro.estimate(season = x[1], nb.af.pride = rep(as.numeric(x[2]), nrow(lions_output_repro)))),
                                                2, median))
repro.nb.af.pride.estimates$lwr = plogis(apply(apply(repro.nb.af.pride.estimates, 
                                                     1,
                                                     FUN = function(x) repro.estimate(season = x[1], nb.af.pride = rep(as.numeric(x[2]), nrow(lions_output_repro)))),
                                               2, FUN = function(x) quantile(x, probs = 0.05)))
repro.nb.af.pride.estimates$upr = plogis(apply(apply(repro.nb.af.pride.estimates,
                                                     1,
                                                     FUN = function(x) repro.estimate(season = x[1], nb.af.pride = rep(as.numeric(x[2]), nrow(lions_output_repro)))),
                                               2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_ReproProb_NbAFPride.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(repro.nb.af.pride.estimates, aes(x = nb.af.pride, y = pred, 
                                        colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = seq(2, 16, 2),
                     labels = seq(2, 16, 2)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  ylim(0, 1) +
  xlab("Number of females\nin the pride") +
  ylab("Reproduction probability\n ") +
  labs(title = "Reproduction probability") +
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


# Average reproduction probability with 
# number of nomadic coalitions in the home range alone

# Empty prediction dataframe
repro.nb.nm.coal.hr.estimates = expand.grid(season = c("wet", "dry"),
                                            nb.nm.coal.hr = (seq(0, 15, 2) - mean(females.data$nb_nm_coal_hr, na.rm = T)) / (2 * sd(females.data$nb_nm_coal_hr, na.rm = T)))

# Fill in prediction and credible intervals
repro.nb.nm.coal.hr.estimates$pred = plogis(apply(apply(repro.nb.nm.coal.hr.estimates, 
                                                        1, 
                                                        FUN = function(x) repro.estimate(season = x[1], nb.nm.coal.hr = rep(as.numeric(x[2]), nrow(lions_output_repro)))),
                                                  2, median))
repro.nb.nm.coal.hr.estimates$lwr = plogis(apply(apply(repro.nb.nm.coal.hr.estimates,
                                                       1,
                                                       FUN = function(x) repro.estimate(season = x[1], nb.nm.coal.hr = rep(as.numeric(x[2]), nrow(lions_output_repro)))), 
                                                 2,
                                                 FUN = function(x) quantile(x, probs = 0.05)))
repro.nb.nm.coal.hr.estimates$upr = plogis(apply(apply(repro.nb.nm.coal.hr.estimates, 
                                                       1,
                                                       FUN = function(x) repro.estimate(season = x[1], nb.nm.coal.hr = rep(as.numeric(x[2]), nrow(lions_output_repro)))),
                                                 2, 
                                                 FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_ReproProb_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(repro.nb.nm.coal.hr.estimates, aes(x = nb.nm.coal.hr, y = pred,
                                          colour = season)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = nb.nm.coal.hr.values,
                     labels = seq(0, 15, 2)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  ylim(0, 1) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Reproduction probability") +
  labs(title = "Reproduction probability") +
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


# Average reproduction probability with number of nomadic coalitions 
# in interaction with the number of females in the pride

# Empty prediction dataframe
repro.nb.nm.coal.hr.nb.af.pride.estimates = expand.grid(season = c("wet", "dry"),
                                                        habitat = c("grassland", "woodland"),
                                                        nb.af.pride = (seq(2, 16, 2) - mean(females.data$nb_af_pride, na.rm = T)) / (2 * sd(females.data$nb_af_pride, na.rm = T)),
                                                        nb.nm.coal.hr = (seq(0, 15, 2) - mean(females.data$nb_nm_coal_hr, na.rm = T)) / (2 * sd(females.data$nb_nm_coal_hr, na.rm = T)))

# Fill in prediction and credible intervals
repro.nb.nm.coal.hr.nb.af.pride.estimates$pred = plogis(apply(apply(repro.nb.nm.coal.hr.nb.af.pride.estimates, 
                                                                    1,
                                                                    FUN = function(x) repro.estimate(season = x[1], 
                                                                                                     habitat = x[2],
                                                                                                     nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_repro)), 
                                                                                                     nb.nm.coal.hr = rep(as.numeric(x[4]),
                                                                                                                           nrow(lions_output_repro)))), 
                                                              2, median))
repro.nb.nm.coal.hr.nb.af.pride.estimates$lwr = plogis(apply(apply(repro.nb.nm.coal.hr.nb.af.pride.estimates, 
                                                                   1, 
                                                                   FUN = function(x) repro.estimate(season = x[1],
                                                                                                    habitat = x[2], 
                                                                                                    nb.af.pride = rep(as.numeric(x[3]),
                                                                                                                      nrow(lions_output_repro)),
                                                                                                    nb.nm.coal.hr = rep(as.numeric(x[4]),
                                                                                                                        nrow(lions_output_repro)))),
                                                             2, 
                                                             FUN = function(x) quantile(x, probs = 0.05)))
repro.nb.nm.coal.hr.nb.af.pride.estimates$upr = plogis(apply(apply(repro.nb.nm.coal.hr.nb.af.pride.estimates,
                                                                   1, 
                                                                   FUN = function(x) repro.estimate(season = x[1],
                                                                                                    habitat = x[2],
                                                                                                    nb.af.pride = rep(as.numeric(x[3]), 
                                                                                                                      nrow(lions_output_repro)), 
                                                                                                    nb.nm.coal.hr = rep(as.numeric(x[4]), 
                                                                                                                        nrow(lions_output_repro)))), 
                                                             2, 
                                                             FUN = function(x) quantile(x, probs = 0.95)))


# Format prediction data for plotting

# Strip labels
temp = c("2 females", "4 females", "6 females", 
         "8 females", "10 females", "12 females")

# Factorize number of females
repro.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride.factor = temp[match(repro.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride, nb.af.pride.values)]
repro.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride.factor = factor(repro.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride.factor, 
                                                                      levels = c("2 females", 
                                                                                 "4 females", 
                                                                                 "6 females", 
                                                                                 "8 females", 
                                                                                 "10 females", 
                                                                                 "12 females"))

# Rename seasons and factorize
repro.nb.nm.coal.hr.nb.af.pride.estimates$season = as.character(repro.nb.nm.coal.hr.nb.af.pride.estimates$season)
repro.nb.nm.coal.hr.nb.af.pride.estimates$season[repro.nb.nm.coal.hr.nb.af.pride.estimates$season == "wet"] = "Wet"
repro.nb.nm.coal.hr.nb.af.pride.estimates$season[repro.nb.nm.coal.hr.nb.af.pride.estimates$season == "dry"] = "Dry"

repro.nb.nm.coal.hr.nb.af.pride.estimates$season = as.factor(repro.nb.nm.coal.hr.nb.af.pride.estimates$season)


# Prediction plot
png(filename = "Predictions_ReproProb_NbNMCoalHRNbAFPride.png",
    width = 6,
    height = 6,
    units = "cm",
    bg = "transparent",
    res = 600,
    type = "cairo")

ggplot(repro.nb.nm.coal.hr.nb.af.pride.estimates[which(repro.nb.nm.coal.hr.nb.af.pride.estimates$habitat == "grassland" & 
                                                       repro.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride %in% 
                                                       c(nb.af.pride.values[1], 
                                                         nb.af.pride.values[4], 
                                                         nb.af.pride.values[6])), ], 
       aes(x = nb.nm.coal.hr, y = pred, color = nb.af.pride.factor, 
           fill = nb.af.pride.factor)) +
  facet_wrap(~ season) +
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  scale_x_continuous(breaks = nb.nm.coal.hr.values,
                     labels = seq(0, 15, 2)) +
  scale_color_manual(name = "Number of females\nin the pride",
                     labels = c("2 females", "8 females", "12 females"),
                     values = viridis::viridis(3, option = "B", end = 0.8)) +
  scale_fill_manual(name = "Number of females\nin the pride",
                    labels = c("2 females", "8 females", "12 females"),
                    values = viridis::viridis(3, option = "B", end = 0.8)) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Reproduction probability") +
  labs(title = "Reproduction probability") + 
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


## 6.2. Recruitment ----
# -----------------

# Subset model output data and remove epsilons and sigmas
lions_output_rec = lions_output_fullGLMM[, grep("rec", 
                                                colnames(lions_output_fullGLMM))]
lions_output_rec.epsilons = lions_output_rec[, grep("epsilon", 
                                                    colnames(lions_output_rec))]
lions_output_rec = lions_output_rec[, - grep("epsilon", 
                                             colnames(lions_output_rec))]
lions_output_rec = lions_output_rec[, - grep("sigma",
                                             colnames(lions_output_rec))]


# Mean, median, and 90% credible interval of the estimates
lions_output_rec_estimates = data.frame(parameter = colnames(lions_output_rec), 
                                        mean = apply(lions_output_rec, 2, mean),
                                        median = apply(lions_output_rec, 2, median),
                                        lowerCI = apply(lions_output_rec, 2, 
                                                        FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_low),
                                        upperCI = apply(lions_output_rec, 2, 
                                                        FUN = function(x) ci(x, method = "ETI", ci = 0.90)$CI_high))

# Full parameter names for plotting
lions_output_rec_estimates$parameter_clean = c("Mean recruitment\n(wet season in the grassland)",
                                               "Mean recruitment\n(dry season in the grassland)",
                                               "Habitat\n(wet season in the grassland)",
                                               "Habitat\n(dry season in the grassland)",
                                               "Habitat\n(wet season in the woodland)",
                                               "Habitat\n(dry seasons in the woodland)",
                                               "Number of females\nin the pride (wet season)",
                                               "Number of females\nin the pride (dry season)",
                                               "Nb females pride:Nb\nnomadic coal HR (wet season)",
                                               "Nb females pride:Nb\nnomadic coal HR (dry season)",
                                               "Number of nomadic coalitions\nin the home range (wet season)",
                                               "Number of nomadic coalitions\nin the home range (dry season)")


# Add season and parameter category
lions_output_rec_estimates$season = rep(c("Wet", "Dry"), 6)
lions_output_rec_estimates$parameter_plot = rep(c("Mean recruitment", 
                                                  "Habitat (grassland)", 
                                                  "Habitat\n(woodland)",
                                                  "Number of females\nin the pride", 
                                                  "Nb females pride:Nb\nnomadic coals HR", 
                                                  "Number of nomadic\ncoalitions in the home range"), 
                                                each = 2)


# Get posterior distribution for density plot
n = nrow(lions_output_fullGLMM) # Number of samples

df.plot.posterior = data.frame(season = rep(rep(c("Wet", "Dry"), each = n), 
                                            length(unique(lions_output_rec_estimates$parameter_plot))),
                               variable = rep(c("Mean recruitment", 
                                                "Habitat (grassland)", 
                                                "Habitat\n(woodland)",
                                                "Number of females\nin the pride",
                                                "Nb females pride:Nb\nnomadic coals HR", 
                                                "Number of nomadic\ncoalitions in the home range"),
                                              each = 2 * n),
                               posterior = c(lions_output_rec[1:nrow(lions_output_rec), ]))

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot.posterior = df.plot.posterior[- which(df.plot.posterior$posterior == 0), ]

# Factorize parameter names and season
df.plot.posterior$variable = factor(df.plot.posterior$variable, 
                                    levels = c("Mean recruitment",
                                               "Number of females\nin the pride",
                                               "Number of nomadic\ncoalitions in the home range",
                                               "Nb females pride:Nb\nnomadic coals HR", "Habitat\n(woodland)")[seq(length(df.plot.posterior$variable), 1, -1)])
df.plot.posterior$season = factor(df.plot.posterior$season,
                                  levels = unique(df.plot.posterior$season))


## 6.2.1. Plotting the effect sizes ----
# ---------------------------------

# Effect size dataframe
df.plot = data.frame(variable = lions_output_rec_estimates$parameter_plot,
                     season = lions_output_rec_estimates$season,
                     mean = lions_output_rec_estimates$mean,
                     median = lions_output_rec_estimates$median,
                     BCI90_upper = lions_output_rec_estimates$upperCI,
                     BCI90_lower = lions_output_rec_estimates$lowerCI,
                     overlapping0 = FALSE)

# Remove 0s to discard habitat effect parameter for the plain, which is the reference level
df.plot = df.plot[- which(df.plot$mean == 0), ]

# Factorize parameter names and season
df.plot$variable = factor(df.plot$variable, 
                          levels = c("Mean recruitment",
                                     "Number of females\nin the pride", 
                                     "Number of nomadic\ncoalitions in the home range", 
                                     "Nb females pride:Nb\nnomadic coals HR", "Habitat\n(woodland)")[length(unique(df.plot$variable)):1])
df.plot$season = factor(df.plot$season, levels = unique(df.plot$season))


# Determine which Bayesian credible intervals overlap 0
for(i in 1:nrow(df.plot)){
  
  if(df.plot$BCI90_upper[i] > 0 & df.plot$BCI90_lower[i] < 0){
    
    df.plot$overlapping0[i] = TRUE
  }
}


# Plot effect size
png(filename = "Recruitment.png", 
    width = 12, 
    height = 10, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

rec_plot = ggplot() +
  ggdist::stat_halfeye(data = df.plot.posterior,  # Density plot
                       aes(x = posterior, y = variable, fill = season),
                       color = NA, size = 1, alpha = 0.4, 
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
  scale_alpha_manual(name = "Significant",
                     labels = c("No", "Yes"),
                     values = c(0.4, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(name = "", breaks = seq(-1, 1, 1), limits = c(-1.4, 1)) +
  theme_general() +
  ggtitle("Recruitment") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

rec_plot

dev.off()


## 6.2.2. Plotting the predictions ----
# --------------------------------

# Function to calculate the predictions
rec.estimate = function(season = "wet",
                        habitat = "grassland",
                        nb.af.pride = 0,
                        nb.nm.coal.hr = 0){
  
  if(season == "wet"){
    
    # Get wet-season parameter values
    mean.rec = lions_output_rec[, grep("mu.rec", 
                                       colnames(lions_output_rec))][, 1]
    beta.nb.af.pride = lions_output_rec[, grep("rec.beta.nb.af.pride", 
                                               colnames(lions_output_rec))][, 1]
    beta.nb.nm.coal.hr = lions_output_rec[, grep("rec.beta.nb.nm.coal.hr", 
                                                 colnames(lions_output_rec))][, 1]
    beta.nb.af.pride.nb.nm.coal.hr = lions_output_rec[, grep("rec.beta.nb.af.pride.nb.nm.coal.hr", 
                                                             colnames(lions_output_rec))][, 1]
    
    # Get wet-season epsilons
    median.epsilons = apply(lions_output_rec.epsilons[, seq(1, 59, 2)], 1, median) 
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_rec))
      
    }
    
    else{
      
      # Get wet-season woodland habitat parameter value
      beta.habitat = lions_output_rec[, grep("rec.beta.habitat.woodland", 
                                             colnames(lions_output_rec))][, 3]
      
    }
  }
  
  else{
    
    # Get dry-season parameter values
    mean.rec = lions_output_rec[, grep("mu.rec", 
                                       colnames(lions_output_rec))][, 2]
    beta.nb.af.pride = lions_output_rec[, grep("rec.beta.nb.af.pride", 
                                               colnames(lions_output_rec))][, 2]
    beta.nb.nm.coal.hr = lions_output_rec[, grep("rec.beta.nb.nm.coal.hr", 
                                                 colnames(lions_output_rec))][, 2]
    beta.nb.af.pride.nb.nm.coal.hr = lions_output_rec[, grep("rec.beta.nb.af.pride.nb.nm.coal.hr", 
                                                             colnames(lions_output_rec))][, 2]
    
    # Get dry-season epsilons
    median.epsilons = apply(lions_output_rec.epsilons[, seq(2, 60, 2)], 1, median)
    
    if(habitat == "grassland"){
      
      # The plain habitat is the reference habitat, the habitat parameter value is thus 0
      beta.habitat = rep(0, nrow(lions_output_rec))
      
    }
    
    else{
      
      # Get dry-season woodland habitat parameter value
      beta.habitat = lions_output_rec[, grep("rec.beta.habitat.woodland", colnames(lions_output_rec))][, 4]
      
    }
  }
  
  # Calculate vital-rate prediction
  pred = mean.rec + beta.habitat + beta.nb.af.pride * nb.af.pride + 
         beta.nb.nm.coal.hr * nb.nm.coal.hr + 
         beta.nb.af.pride.nb.nm.coal.hr * nb.af.pride * nb.nm.coal.hr + 
         median.epsilons
  
  return(pred)
  
}


# Average survival per season and habitat

# Empty prediction dataframe
rec.mean.season.estimates = expand.grid(season = c("wet", "dry"),
                                        habitat = c("grassland", "woodland"))

# Fill in prediction and credible intervals
rec.season.estimates$pred = exp(apply(apply(rec.mean.season.estimates, 
                                            1, 
                                            FUN = function(x) rec.estimate(season = x[1], habitat = x[2])), 
                                      2, median))
rec.season.estimates$lwr = exp(apply(apply(rec.mean.season.estimates,
                                           1, 
                                           FUN = function(x) rec.estimate(season = x[1], habitat = x[2])), 
                                     2, FUN = function(x) quantile(x, probs = 0.05)))
rec.season.estimates$upr = exp(apply(apply(rec.mean.season.estimates,
                                           1, 
                                           FUN = function(x) rec.estimate(season = x[1], habitat = x[2])), 
                                     2, FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_Recruit_SeasonHabitat.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rec.season.estimates, aes(x = season, y = pred, colour = habitat)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.8), width = 0.5) +
  scale_x_discrete(labels = c("Wet", "Dry")) +
  scale_color_manual(name = "Habitat",
                     labels = c("Grassland", "Woodland"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Season") +
  ylab("Recruitment") +
  labs(title = "Recruitment") + 
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


# Average recruitment with number of females in the pride

# Empty prediction dataframe
rec.nb.af.pride.estimates = expand.grid(season = c("wet", "dry"), 
                                        nb.af.pride = (seq(1, 16, 1) - mean(females.data$nb_af_pride, na.rm = T)) / (2 * sd(females.data$nb_af_pride, na.rm = T)))

# Fill in prediction and credible intervals
rec.nb.af.pride.estimates$pred = exp(apply(apply(rec.nb.af.pride.estimates, 
                                                 1,
                                                 FUN = function(x) rec.estimate(season = x[1], 
                                                                                nb.af.pride = rep(as.numeric(x[2]), nrow(lions_output_rec)))),
                                           2, median))
rec.nb.af.pride.estimates$lwr = exp(apply(apply(rec.nb.af.pride.estimates,
                                                1,
                                                FUN = function(x) rec.estimate(season = x[1],
                                                                               nb.af.pride = rep(as.numeric(x[2]), nrow(lions_output_rec)))),
                                          2,
                                          FUN = function(x) quantile(x, probs = 0.05)))
rec.nb.af.pride.estimates$upr = exp(apply(apply(rec.nb.af.pride.estimates,
                                                1, 
                                                FUN = function(x) rec.estimate(season = x[1], 
                                                                               nb.af.pride = rep(as.numeric(x[2]), nrow(lions_output_rec)))), 
                                          2, 
                                          FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Models_Output/Plots/Predictions_Recruit_SeasonalCoeffs_NbAFPride.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rec.nb.af.pride.estimates, aes(x = nb.af.pride, y = pred, 
                                      colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = nb.af.pride.values,
                     labels = seq(1, 16, 1)) +
  scale_color_manual(name = "Season",
                     labels = c("Dry", "Wet"),
                     values = c(cbbPalette[2], cbbPalette[4])) +
  scale_fill_manual(name = "Season",
                    labels = c("Dry", "Wet"),
                    values = c(cbbPalette[2], cbbPalette[4])) +
  xlab("Number of adult females\nin the pride") +
  ylab("Recruitment") +
  labs(title = "Recruitment") + 
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


# Average recruitment with number of females in the pride
# in interaction with the number of nomadic coalitions in the home range

# Empty prediction dataframe 
rec.nb.af.pride.nb.nm.coal.hr.estimates = expand.grid(season = c("wet", "dry"),
                                                      habitat = c("grassland", "woodland"),
                                                      nb.af.pride = (seq(2, 12, 2) - mean(females.data$nb_af_pride, na.rm = T)) / (2 * sd(females.data$nb_af_pride, na.rm = T)),
                                                      nb.nm.coal.hr = (seq(0, 15, 2) - mean(females.data$nb_nm_coal_hr, na.rm = T)) / (2 * sd(females.data$nb_nm_coal_hr, na.rm = T)))

# Fill in prediction and credible intervals
rec.nb.af.pride.nb.nm.coal.hr.estimates$pred = exp(apply(apply(rec.nb.af.pride.nb.nm.coal.hr.estimates, 
                                                               1, 
                                                               FUN = function(x) rec.estimate(season = x[1], 
                                                                                              habitat = x[2], 
                                                                                              nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_rec)), 
                                                                                              nb.nm.coal.hr = rep(as.numeric(x[4]), nrow(lions_output_rec)))), 
                                                         2, median))
rec.nb.af.pride.nb.nm.coal.hr.estimates$lwr = exp(apply(apply(rec.nb.af.pride.nb.nm.coal.hr.estimates, 
                                                              1,
                                                              FUN = function(x) rec.estimate(season = x[1],
                                                                                             habitat = x[2], 
                                                                                             nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_rec)),
                                                                                             nb.nm.coal.hr = rep(as.numeric(x[4]), nrow(lions_output_rec)))), 
                                                        2, FUN = function(x) quantile(x, probs = 0.05)))
rec.nb.af.pride.nb.nm.coal.hr.estimates$upr = exp(apply(apply(rec.nb.af.pride.nb.nm.coal.hr.estimates,
                                                              1, 
                                                              FUN = function(x) rec.estimate(season = x[1],
                                                                                             habitat = x[2], 
                                                                                             nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_rec)),
                                                                                             nb.nm.coal.hr = rep(as.numeric(x[4]), nrow(lions_output_rec)))),
                                                        2, FUN = function(x) quantile(x, probs = 0.95)))


# Format prediction data for plotting

# Strip labels
temp = c("No nomadic coalition", "2 nomadic coalitions", "4 nomadic coalitions", 
         "6 nomadic coalitions", "8 nomadic coalitions", "10 nomadic coalitions", 
         "12 nomadic coalitions", "14 nomadic coalitions")

# Factorize numbeer of nomadic coalitions
rec.nb.af.pride.nb.nm.coal.hr.estimates$nb.nm.coal.hr.factor = temp[match(rec.nb.af.pride.nb.nm.coal.hr.estimates$nb.nm.coal.hr, 
                                                                          nb.nm.coal.hr.values)]
rec.nb.af.pride.nb.nm.coal.hr.estimates$nb.nm.coal.hr.factor = factor(rec.nb.af.pride.nb.nm.coal.hr.estimates$nb.nm.coal.hr.factor, 
                                                                      levels = c("No nomadic coalition", 
                                                                                 "2 nomadic coalitions", 
                                                                                 "4 nomadic coalitions", 
                                                                                 "6 nomadic coalitions",
                                                                                 "8 nomadic coalitions", 
                                                                                 "10 nomadic coalitions", 
                                                                                 "12 nomadic coalitions", 
                                                                                 "14 nomadic coalitions"))


# Rename seasons and factorize
rec.nb.af.pride.nb.nm.coal.hr.estimates$season = as.character(rec.nb.af.pride.nb.nm.coal.hr.estimates$season)
rec.nb.af.pride.nb.nm.coal.hr.estimates$season[rec.nb.af.pride.nb.nm.coal.hr.estimates$season == "wet"] = "Wet"
rec.nb.af.pride.nb.nm.coal.hr.estimates$season[rec.nb.af.pride.nb.nm.coal.hr.estimates$season == "dry"] = "Dry"

rec.nb.af.pride.nb.nm.coal.hr.estimates$season = as.factor(rec.nb.af.pride.nb.nm.coal.hr.estimates$season)


# Prediction plot
png(filename = "Predictions_Recruit_NbAFPrideNbNMCoalHR.png",
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rec.nb.af.pride.nb.nm.coal.hr.estimates[which(rec.nb.af.pride.nb.nm.coal.hr.estimates$habitat == "grassland" & 
                                                     rec.nb.af.pride.nb.nm.coal.hr.estimates$nb.nm.coal.hr %in% 
                                                     c(nb.nm.coal.hr.values[1], 
                                                       nb.nm.coal.hr.values[4],
                                                       nb.nm.coal.hr.values[6],
                                                       nb.nm.coal.hr.values[8])), ], 
       aes(x = nb.af.pride, y = pred, color = nb.nm.coal.hr.factor, 
           fill = nb.nm.coal.hr.factor)) +
  facet_wrap(~ season) +
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  scale_x_continuous(breaks = nb.af.pride.values,
                     labels = seq(2, 12, 2)) +
  scale_color_manual(name = "Number of nomadic coalitions\nin the home range",
                     labels = c("No coalition", "6 coalitions", 
                                "10 coalitions", "14 coalitions"),
                     values = viridis::viridis(4, option = "B", end = 0.8)) +
  scale_fill_manual(name = "Number of nomadic coalitions\nin the home range",
                    labels = c("No coalition", "6 coalitions", 
                               "10 coalitions", "14 coalitions"),
                    values = viridis::viridis(4, option = "B", end = 0.8)) +
  xlab("Number of adut females\nin the pride") +
  ylab("Recruitment") +
  labs(title = "Recruitment") + 
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


# Average recruitment with number of nomadic coalitions in the home range

# Empty prediction dataframe
rec.nb.nm.coal.hr.estimates = expand.grid(season = c("wet", "dry"), 
                                          nb.nm.coal.hr = (seq(0, 15, 2) - mean(females.data$nb_nm_coal_hr, na.rm = T)) / (2 * sd(females.data$nb_nm_coal_hr, na.rm = T)))

# Fill in prediction and credible intervals
rec.nb.nm.coal.hr.estimates$pred = exp(apply(apply(rec.nb.nm.coal.hr.estimates, 
                                                   1, 
                                                   FUN = function(x) rec.estimate(season = x[1], 
                                                                                  nb.nm.coal.hr = rep(as.numeric(x[2]), nrow(lions_output_rec)))),
                                             2, mean))
rec.nb.nm.coal.hr.estimates$lwr = exp(apply(apply(rec.nb.nm.coal.hr.estimates,
                                                  1, 
                                                  FUN = function(x) rec.estimate(season = x[1], 
                                                                                 nb.nm.coal.hr = rep(as.numeric(x[2]), nrow(lions_output_rec)))), 
                                            2, 
                                            FUN = function(x) quantile(x, probs = 0.05)))
rec.nb.nm.coal.hr.estimates$upr = exp(apply(apply(rec.nb.nm.coal.hr.estimates,
                                                  1,
                                                  FUN = function(x) rec.estimate(season = x[1], 
                                                                                 nb.nm.coal.hr = rep(as.numeric(x[2]), nrow(lions_output_rec)))),
                                            2,
                                            FUN = function(x) quantile(x, probs = 0.95)))


# Prediction plot
png(filename = "Predictions_Recruit_NbNMCoalHR.png", 
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rec.nb.nm.coal.hr.estimates, aes(x = nb.nm.coal.hr, y = pred, 
                                        colour = season)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) +
  scale_x_continuous(breaks = nb.nm.coal.hr.values,
                     labels = seq(0, 15, 2)) +
  scale_color_manual(name = "Season",
                     labels = c("Wet", "Dry"),
                     values = c(cbbPalette[4], cbbPalette[2])) +
  scale_fill_manual(name = "Season",
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[4], cbbPalette[2])) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Recruitment\n ") +
  labs(title = "Recruitment") + 
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


# Average recruitment with number of nomadic coalitions
# in interaction with the number of females in the pride

# Empty prediction dataframe
rec.nb.nm.coal.hr.nb.af.pride.estimates = expand.grid(season = c("wet", "dry"),
                                                      habitat = c("grassland", "woodland"),
                                                      nb.af.pride = (seq(2, 12, 2) - mean(females.data$nb_af_pride, na.rm = T)) / (2 * sd(females.data$nb_af_pride, na.rm = T)),
                                                      nb.nm.coal.hr = (seq(0, 15, 2) - mean(females.data$nb_nm_coal_hr, na.rm = T)) / (2 * sd(females.data$nb_nm_coal_hr, na.rm = T)))

# Fill in prediction and credible intervals
rec.nb.nm.coal.hr.nb.af.pride.estimates$pred = exp(apply(apply(rec.nb.nm.coal.hr.nb.af.pride.estimates,
                                                               1, 
                                                               FUN = function(x) rec.estimate(season = x[1], 
                                                                                              habitat = x[2],
                                                                                              nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_rec)),
                                                                                              nb.nm.coal.hr = rep(as.numeric(x[4]), nrow(lions_output_rec)))),
                                                         2, median))
rec.nb.nm.coal.hr.nb.af.pride.estimates$lwr = exp(apply(apply(rec.nb.nm.coal.hr.nb.af.pride.estimates, 
                                                              1, 
                                                              FUN = function(x) rec.estimate(season = x[1],
                                                                                             habitat = x[2], 
                                                                                             nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_rec)), 
                                                                                             nb.nm.coal.hr = rep(as.numeric(x[4]), nrow(lions_output_rec)))), 
                                                        2, FUN = function(x) quantile(x, probs = 0.05)))
rec.nb.nm.coal.hr.nb.af.pride.estimates$upr = exp(apply(apply(rec.nb.nm.coal.hr.nb.af.pride.estimates, 
                                                              1, 
                                                              FUN = function(x) rec.estimate(season = x[1], 
                                                                                             habitat = x[2],
                                                                                             nb.af.pride = rep(as.numeric(x[3]), nrow(lions_output_rec)), 
                                                                                             nb.nm.coal.hr = rep(as.numeric(x[4]), nrow(lions_output_rec)))), 
                                                        2, 
                                                        FUN = function(x) quantile(x, probs = 0.95)))


# Format prediction data for plotting

# Strip labels
temp = c("2 females", "4 females", "6 females",
         "8 females", "10 females", "12 females")

# Factorize number of females
rec.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride.factor = temp[match(rec.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride, 
                                                                        nb.af.pride.values)]
rec.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride.factor = factor(rec.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride.factor, 
                                                                    levels = c("2 females", 
                                                                               "4 females", 
                                                                               "6 females", 
                                                                               "8 females", 
                                                                               "10 females",
                                                                               "12 females"))

# Rename seasons and factorize
rec.nb.nm.coal.hr.nb.af.pride.estimates$season = as.character(rec.nb.nm.coal.hr.nb.af.pride.estimates$season)
rec.nb.nm.coal.hr.nb.af.pride.estimates$season[rec.nb.nm.coal.hr.nb.af.pride.estimates$season == "wet"] = "Wet"
rec.nb.nm.coal.hr.nb.af.pride.estimates$season[rec.nb.nm.coal.hr.nb.af.pride.estimates$season == "dry"] = "Dry"

rec.nb.nm.coal.hr.nb.af.pride.estimates$season = as.factor(rec.nb.nm.coal.hr.nb.af.pride.estimates$season)


# Prediction plot
png(filename = "Predictions_Recruit_NbNMCoalHRNbAFPride.png",
    width = 6, 
    height = 6, 
    units = "cm", 
    bg = "transparent", 
    res = 600, 
    type = "cairo")

ggplot(rec.nb.nm.coal.hr.nb.af.pride.estimates[which(rec.nb.nm.coal.hr.nb.af.pride.estimates$habitat == "grassland" & 
                                                     rec.nb.nm.coal.hr.nb.af.pride.estimates$nb.af.pride %in% 
                                                     c(nb.af.pride.values[1], 
                                                       nb.af.pride.values[4], 
                                                       nb.af.pride.values[6])), ], 
       aes(x = nb.nm.coal.hr, y = pred, color = nb.af.pride.factor, 
           fill = nb.af.pride.factor)) +
  facet_wrap(~ season) +
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  scale_x_continuous(breaks = nb.nm.coal.hr.values,
                     labels = seq(0, 15, 2)) +
  scale_color_manual(name = "Number of females\nin the pride",
                     labels = c("2 females", "8 females", "12 females"),
                     values = viridis::viridis(3, option = "B", end = 0.8)) +
  scale_fill_manual(name = "Number of females\nin the pride",
                    labels = c("2 females", "8 females", "12 females"),
                    values = viridis::viridis(3, option = "B", end = 0.8)) +
  xlab("Number of nomadic coalitions\nin the home range") +
  ylab("Recruitment") +
  labs(title = "Recruitment") + 
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








###########################################################################
#
# 7. Epsilon posterior distributions ----
#
###########################################################################

# Dataframe to match full parameter names, model output column, and file name
# We add empty slots for plot aesthetics reasons.

params.labels = data.frame(# Output model column name
                           col.label = c("repro",
                                         "rec", 
                                         "1", 
                                         "2", 
                                         "3", 
                                         "4", 
                                         "5", 
                                         "6", 
                                         "7", 
                                         "8", 
                                         "9", 
                                         "10", 
                                         "11", 
                                         "12", 
                                         "13"),
                           # Full parameter name
                           param.name = c("Reproduction probability",
                                          "Recruitment",
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
                           # Fil name
                           param.file = c("ReproProb",
                                          "Recruitment", 
                                          "1", 
                                          "2", 
                                          "3", 
                                          "4", 
                                          "5", 
                                          "6", 
                                          "7", 
                                          "8", 
                                          "9", 
                                          "10", 
                                          "11", 
                                          "12", 
                                          "13"))

# Epsilon distribution dataframe
nSamples = 1

epsilon_df = data.frame(parameter = rep(rep(params.labels$param.name, 
                                            each = nSamples * length(seq(1985, 2014))), 2),
                        label = rep(rep(params.labels$col.label, 
                                        each = nSamples * length(seq(1985, 2014))), 2),
                        season = rep(c("Wet", "Dry"),
                                     each = nSamples * length(seq(1985, 2014)) * length(params.labels$param.name)),
                        year = rep(rep(rep(seq(1985, 2014),
                                           each = nSamples), 
                                       length(params.labels$param.name)), 2), 
                        median = NA, 
                        lwr = NA,
                        upr = NA)


# Fill in epsilon distribution dataframe
epsilon_df$mean[which(is.na(epsilon_df$mean))] = rnorm(length(epsilon_df$mean[which(is.na(epsilon_df$mean))]), 0, 4)
epsilon_df$lwr[which(is.na(epsilon_df$lwr))] = rnorm(length(epsilon_df$lwr[which(is.na(epsilon_df$lwr))]), 0, 4)
epsilon_df$upr[which(is.na(epsilon_df$upr))] = rnorm(length(epsilon_df$upr[which(is.na(epsilon_df$upr))]), 0, 4)

for(p in 1:nrow(params.labels[1:2, ])){
  
  param = params.labels[p, ] # Get parameter full name
  
  epsilon.colname = paste0("epsilon.", param[1]) # Get model output column name
  
  # Subset wet-season epsilons
  epsilons_wet = lions_output_fullGLMM[, grep(epsilon.colname,
                                              colnames(lions_output_fullGLMM))] 
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
  
  # Subset dry-season epsilons
  epsilons_dry = lions_output_fullGLMM[, grep(epsilon.colname, 
                                              colnames(lions_output_fullGLMM))] 
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
  ggplot(epsilon_df, aes(x = year, y = mean, colour = season)) +
  facet_wrap(~ parameter, ncol = 3, scales = "free_y") +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), alpha = 0.2) + 
  scale_fill_manual(name = "Season",
                    values = c(cbbPalette[4], cbbPalette[2]),
                    labels = c("Wet", "Dry")) +
  scale_colour_manual(name = "Season",
                      values = c(cbbPalette[4], cbbPalette[2]),
                      labels = c("Wet", "Dry")) +
  xlab("Year") +
  ylab("Epsilon") + 
  theme_minimal() %+replace%    
  theme(axis.text = element_text(colour = colPlot, size = 7, family = font),
        axis.title.x = element_text(colour = colPlot, size = 8, family = font, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(colour = colPlot, size = 8, family = font, margin = margin(t = 0, r = 15, b = 0, l = 0), angle = 90),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = colPlot),
        axis.ticks = element_line(colour = colPlot),
        plot.background = element_rect(fill = colBG, color = colBG, size = 0),
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
