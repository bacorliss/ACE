

# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)

p_load(ggplot2)
p_load(tibble)
p_load(broom)
p_load(tidyr)
p_load(dplyr)
p_load(boot)
p_load(readr)
p_load(gplots)
# User defined libraries
source("R/mmd.R")
source("R/equivalence_contests.R")


# Figure parameters
fig_num = "9" 
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/F",fig_num, "/",sep="")
# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_sims = 1e3
n_samples = 1e3
n_obs = 50
rand.seed = 1

parallel_sims = TRUE
include_bf = TRUE
# scale_contest_path = paste("figure/F", fig_num, "/F", fig_num,"_scale_contest_results.csv",sep="")

# Test how each metric responds to changing each characteristic of similarity





# Pearson Correlation between abs(mu_d1) versus mean value of stats
# sigmad = 1
#------------------------------------------------------------------------------
set.seed(rand.seed)
mu_d_vect = seq(-10,10,0.5)
indvar_vect = 1+mu_d_vect; n_sims = length(indvar_vect)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1a_stat_correlation_mud_sweep", sep = "")
df_init_1a <- generateExperiment_Data(n_samples, n_obs, n_sims = n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = sqrt((1^2)/2), 
                                   mus_1b  = indvar_vect, 
                                   sigmas_1b = sqrt((1^2)/2),
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(0,n_sims),  
                                   sigmas_2b = 1, 
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize_1a <- process_esize_simulations(df_init_1a, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize_1a$df_es, indvar = "mu_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])




# Pearson Correlation between abs(mu_d1) versus mean value of stats
# sigmad = 1
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigma_d_vect = seq(.1,10,0.1); n_sims = length(sigma_d_vect)
sigma_ab_vect = sqrt((sigma_d_vect^2)/2)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1b_stat_correlation_sigmad_sweep_mu-0", sep = "")
df_init_1a <- generateExperiment_Data(n_samples, n_obs, n_sims = n_sims, rand.seed, 
                                      mus_1a  = rep(0,n_sims), 
                                      sigmas_1a = sigma_ab_vect, 
                                      mus_1b  = rep(0,n_sims), 
                                      sigmas_1b = sigma_ab_vect,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = sqrt((1^2)/2), 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = sqrt((1^2)/2), 
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize_1a <- process_esize_simulations(df_init_1a, gt_colname = gt_colnames, 
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize_1a$df_es, indvar = "sigma_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])



# Pearson Correlation between abs(mu_d1) versus mean value of stats
# sigmad = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigma_d_vect = seq(.1,10,0.1); n_sims = length(sigma_d_vect)
sigma_ab_vect = sqrt((sigma_d_vect^2)/2)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1b_stat_correlation_sigmad_sweep_mu-10", sep = "")
df_init_1a <- generateExperiment_Data(n_samples, n_obs, n_sims = n_sims, rand.seed, 
                                      mus_1a  = rep(0,n_sims), 
                                      sigmas_1a = sigma_ab_vect, 
                                      mus_1b  = rep(10,n_sims), 
                                      sigmas_1b = sigma_ab_vect,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = sqrt((1^2)/2), 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = sqrt((1^2)/2), 
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize_1a <- process_esize_simulations(df_init_1a, gt_colname = gt_colnames, 
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize_1a$df_es, indvar = "sigma_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])




# Change sigma at mu=10, calculate pearson for each stat across sigma




# Change rmu w/ sigma=1, calculate pearson for each stat across rmu


# Change rsigma at mu=0, calculate pearson for each stat across rsigma


# Change rsigma at mu=10, calculate pearson for each stat across rsigma


# Summarizing table showing MMD is the only stat with a solid color that is statistically significant from







# CHange sampling with fixed std and mu