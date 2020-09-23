

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




################################################################################

# Unscaled Space________________________________________________________________




# Mu
# Pearson rho of abs(mu_d1) versus abs(mean of each stat)
# Sweep mu_1d, sigma_d = 1
#------------------------------------------------------------------------------
# Fixed mu_d, but as it increases, rmu_d decreases
set.seed(rand.seed)
mus_d_vect = seq(0.1,5,0.1)
mus_a_vect = mus_d_vect*2
mus_b_vect = mus_d_vect + mus_a_vect
rmus_d_vect = mus_d_vect/mus_a_vectt

n_sims = length(mus_b_vect)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1a_stat_correlation_mud_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples = n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_a_vect, 
                                   sigmas_1a = sqrt((1^2)/2), 
                                   mus_1b  = mus_b_vect, 
                                   sigmas_1b = sqrt((1^2)/2),
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(0,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)

df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "mu_1md", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])


# Sigma
# Pearson rho of sigma versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigmas_d_vect = seq(0.1,  10,  0.1); n_sims = length(sigmas_d_vect)
mus_ab_vect =   seq(0.1,  10,  0.1)
rsigmas_d_vect = sigmas_d_vect/mus_ab_vect
sigmas_ab_vect = sqrt((sigmas_d_vect^2)/2)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1b_stat_correlation_sigmad_sweep_mu-10", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_ab_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_ab_vect, 
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "sigma[DM]",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "sigma_1md", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])






# Relative Mean
# rmu_d1 vs. mean of stats, pearson rho
# sigmad = 1
#------------------------------------------------------------------------------
# Fixed mu_d, but as it increases, rmu_d decreases
source("R/equivalence_contests.R")
set.seed(rand.seed)
rmus_d_vect = seq(0.25, 4,0.1); n_sims = length(rmus_d_vect)
mus_a_vect = 1; 
mus_b_vect = rmus_d_vect*mus_a_vect
mus_d_vect = mus_b_vect -mus_a_vect



sigmas_d_vect = mus_b_vect
rsigmas_d_vect = mus_b_vect/mus_a_vect
  
sigmas_ab_vect = sqrt(0.5 * sigmas_d_vect^2)


gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1a_stat_correlation_mud_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed,
                                      mus_1a  = mus_a_vect,
                                      sigmas_1a = sigmas_ab_vect,
                                      mus_1b  = mus_b_vect,
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(0,n_sims),
                                      sigmas_2a = 1,
                                      mus_2b  = rep(0,n_sims),
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames,
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <-
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "rmu_1md",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])

# Pearson rho of abs(mu_d1) versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigmas_d_vect = seq(0.1,  10,  0.1); n_sims = length(sigmas_d_vect)
mus_ab_vect =   seq(0.1,  10,  0.1)
rsigmas_d_vect = sigmas_d_vect/mus_ab_vect
sigmas_ab_vect = sqrt((sigmas_d_vect^2)/2)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1b_stat_correlation_sigmad_sweep_mu-10", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_ab_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_ab_vect, 
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "sigma_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])


# rigma_d versus abs(stat), Pearson Rho
# mu_1a = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigmas_d_vect = seq(0.1,  10,  0.1); n_sims = length(sigmas_d_vect)
mus_ab_vect =   seq(0.1,  10,  0.1)
rsigmas_d_vect = sigmas_d_vect/mus_ab_vect
sigmas_ab_vect = sqrt((sigmas_d_vect^2)/2)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1b_stat_correlation_sigmad_sweep_mu-10", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_ab_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_ab_vect, 
                                      sigmas_1b = sigmas_ab_vect+10,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "sigma_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])




# Sample Size:   Pearson rho of sigma versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
n_1ab_vect = seq(100, 16,  -2); n_sims = length(n_1ab_vect)
mus_ab_vect = 1
sigmas_ab_vect = 1

fig_name = paste("F", fig_num, "_1b_stat_correlation_n_sweep_mu-1", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_ab_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_ab_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(0,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_1ab_vect, n_1b = n_1ab_vect, n_2a = 30, n_2b = 30,
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "sigma[DM]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
source("R/equivalence_contests.R")
df_mud_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "df_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,
                           stats_basenames = effect_size_dict[[2]],
                           stats_labels = effect_size_dict[[4]])


