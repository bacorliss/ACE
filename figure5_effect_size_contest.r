

### Figure 2: MMD outperforms previous statistics of effect size.
#
# Metrics of sample effect size
#
##############################################
#  Un-standardized measures of effect size
#   * difference of sample means
#   * relative difference of means
#   * difference of sample variances
#  Standardized measures of effect size
#   * cohen's d
#   * hedges g
#   * glasses delta
#   * MHD

library(pacman)
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(broom)
p_load(tidyr)
p_load(dplyr)
p_load(boot)
# User defined libraries
source("R/mmd.R")
source("R/effect_size_contests.R")


boot_xbar <- function(x, ind)  mean(x[ind])
extend_max_lim <- function(x,prc) max(x) + prc* (max(x)-min(x))

# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_sims = 1e3
n_samples = 1e2
n_obs = 50
rand.seed = 0
fig_basename = "F5"
parallel_sims = TRUE
include_bf = TRUE
switch_group_id = TRUE
switch_mean_sign = TRUE

###############################################################################
#
# Un-transformed Error

# Contest 1-1) Quantify Error rate with each metric predicting experiment with
# Lower mean difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_1a_esize_", "contest_mu_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,-1,1), 
                                   sigmas_1d = rep(0.02,n_sims),
                                   mus_2a  = rep(10,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,-1,1), 
                                   sigmas_2d = rep(0.02,n_sims),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = switch_mean_sign,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                     y_ax_str = "abs(~mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))
# [ Far from zero ]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_1b_esize_contest_mu_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,5,10), 
                                    sigmas_1d = rep(0.01,n_sims),
                                   mus_2a  = rep(10,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,5,10), 
                                    sigmas_2d = rep(0.01,n_sims),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = switch_mean_sign,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                 y_ax_str = "abs(~mu[D]*phantom(.))",
                                 include_bf = include_bf, parallel_sims = parallel_sims,
                                 fig_name = paste(fig_name, ".tiff",sep = ""))




# COntest 2) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_2a_esize_", "contest_sigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(0.1,n_sims), 
                                   sigmas_1d = runif(n_sims, 1, 20),
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(0.1,n_sims), 
                                   sigmas_2d = runif(n_sims, 1, 20),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = switch_mean_sign,
                                   fig_name = paste(fig_name, ".tiff",sep = "")) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))
# [Far from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_2b_esize_contest_sigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(10,n_sims), 
                                   sigmas_1d = runif(n_sims, 5, 20),
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(10,n_sims), 
                                   sigmas_2d = runif(n_sims, 5, 20),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = switch_mean_sign,
                                   fig_name = paste(fig_name, ".tiff",sep = "")) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))



# Contest 3) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means and mean difference in means [Free form]
# [Near from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_3a_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,1,1), 
                                   sigmas_1a = runif(n_sims,1,1), 
                                   mus_1d  = runif(n_sims,-1,1), 
                                   sigmas_1d = runif(n_sims,1,5),
                                   mus_2a  = runif(n_sims,1,1),  
                                   sigmas_2a = runif(n_sims,1,1),
                                   mus_2d  = runif(n_sims,-1,1), 
                                   sigmas_2d = runif(n_sims,1,5),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = switch_mean_sign,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                     y_ax_str = "abs(~mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_mu.tiff",sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_sigma.tiff",sep = ""))
# [Far from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_3b_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims, 1, 1), 
                                   sigmas_1a = runif(n_sims, 1, 1), 
                                   mus_1d  = runif(n_sims, 5, 10), 
                                   sigmas_1d = runif(n_sims, 5, 50),
                                   mus_2a  = runif(n_sims, 1, 1),  
                                   sigmas_2a = runif(n_sims, 1, 1),
                                   mus_2d  = runif(n_sims, 5, 10), 
                                   sigmas_2d = runif(n_sims, 5, 50),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = switch_mean_sign,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                     y_ax_str = "abs(~mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_mu.tiff",sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_sigma.tiff",sep = ""))



###############################################################################
#
# Realative Error
#
##############################################################################

# Contest 4) Quantify Error rate with each metric predicting experiment with
# Lower relative mean difference in means
# [Near zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_4a_esize_contest_rmu_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,-5,10), 
                                   sigmas_1d = rep(0.01,n_sims),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,-25, 50), 
                                   sigmas_2d = rep(0.05,n_sims),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "abs(~r*mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))
# [ Far from zero ]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_4b_esize_contest_rmu_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,10,40), 
                                   sigmas_1d = rep(0.2,n_sims),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,50, 200), 
                                   sigmas_2d = rep(1,n_sims),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "abs(~r*mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))


# Contest 5) Quantify Error rate with each metric predicting experiment with
# Lower relative STD of difference in means
# [u_d near zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_5a_esize_contest_rsigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1, n_sims),
                                   mus_1d  = runif(n_sims, -1, 1), 
                                   sigmas_1d = runif(n_sims, 1, 20),
                                   mus_2a  = rep(40,n_sims), 
                                   sigmas_2a = rep(1, n_sims),
                                   mus_2d  =  runif(n_sims, -1, 1), 
                                   sigmas_2d = runif(n_sims, 4, 80),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = "")) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                                     y_ax_str = "r*sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))
# [ Far from zero ]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_5b_esize_contest_rsigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(10,n_sims), 
                                    sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(20,n_sims), 
                                    sigmas_1d = runif(n_sims, 10, 50),
                                   mus_2a  = rep(40,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(80,n_sims), 
                                    sigmas_2d = runif(n_sims, 40, 200),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = "")) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                      y_ax_str = "r*sigma[d]",
                      include_bf = include_bf, parallel_sims = parallel_sims,
                      fig_name = paste(fig_name, ".tiff",sep = ""))





# Contest 6) Quantify Error rate with each metric predicting experiment with
# Lower rel. STD, and rel. mean, of difference in means [Free form]
# [Near from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_6a_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,1,1), 
                                   sigmas_1a = runif(n_sims,1,1), 
                                   mus_1d  = runif(n_sims,-1,1), 
                                   sigmas_1d = runif(n_sims,1,5),
                                   mus_2a  = runif(n_sims,1,1),  
                                   sigmas_2a = runif(n_sims,1,1),
                                   mus_2d  = runif(n_sims,-1,1), 
                                   sigmas_2d = runif(n_sims,1,5),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""))  
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "abs(~r*mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rmu.tiff",sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "r*sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rsigma.tiff",sep = ""))
# [Far from zero]
#
#------------------------------------------------------------------------------
fig_name = paste(fig_basename, "_6b_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims, 10 ,10), 
                                   sigmas_1a = runif(n_sims, 1, 1), 
                                   mus_1d  = runif(n_sims, 5, 15), 
                                   sigmas_1d = runif(n_sims, 15, 60),
                                   
                                   mus_2a  = runif(n_sims, 20, 20),  
                                   sigmas_2a = runif(n_sims, 1, 1),
                                   mus_2d  = runif(n_sims, 10, 30), 
                                   sigmas_2d = runif(n_sims, 30, 120),
                                   switch_group_id = switch_group_id,
                                   switch_mean_sign = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = "")) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "abs(~r*mu[D]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rmu.tiff",sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                                     y_ax_str = "r*sigma[d]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rsigma.tiff",sep = ""))












# Contest 7) which metric is best at discerning experiments below a threshold 
# difference in means
#
#------------------------------------------------------------------------------




# Contest 8) which metric is best at discerning experiments below a threshold 
# relative difference in means
#
#------------------------------------------------------------------------------




