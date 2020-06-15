

### Figure 2: MMD outperforms previous statistics of effect size.
#
# Metrics of sample effect size
#
##############################################
#  Unstandardized measures of effect size
#   * difference of sample means
#   * relative difference of means
#   * difference of sample variances
#  Standardized measures of effect size
#   * cohen's d
#   * hedges g
#   * glasses delta
#   * MHD


library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
# library(gridExtra)
# library(grid)
library(tidyr)
library(cowplot)
library(dplyr)
# library(effsize)
library(boot)
# User defined libraries
source("R/mmd.R")
source("R/effect_size_contests.R")


boot_xbar <- function(x, ind)  mean(x[ind])
extend_max_lim <- function(x,prc) max(x) + prc* (max(x)-min(x))

# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomnly chosen
n_sims = 1e2
n_samples = 1e3
n_obs = 50
rand.seed = 0
fig_basename = "f_5"





###############################################################################
#
# Untrasformed Error

# Contest 1-1) Quantify Error rate with each metric predicting experiment with
# Lower mean difference in means
# [ Far from zero ]
#

#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,5,10), 
                                    sigmas_1d = rep(0,n_sims),
                                   mus_2a  = rep(1,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,5,10), 
                                    sigmas_2d = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                 y_ax_str = "mu[d]",
                                 fig_name = paste(fig_basename, "_1a_esize_",
                                                  "contest_mu_far_zero.tiff",
                                                  sep = ""))
# [Near from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,-1,1), 
                                    sigmas_1d = rep(0,n_sims),
                                   mus_2a  = rep(1,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,-1,1), 
                                    sigmas_2d = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                     y_ax_str = "mu[d]",
                                     fig_name = paste(fig_basename, "_1b_esize_",
                                                      "contest_mu_near_zero.tiff", 
                                                      sep = ""))


# COntest 2) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means
# [Far from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(99,n_sims), 
                                   sigmas_1d = runif(n_sims, 1, 20),
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(99,n_sims), 
                                   sigmas_2d = runif(n_sims, 1, 20))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     fig_name = paste(fig_basename, "_2a_esize_",
                                                      "contest_sigma_far_zero.tiff", 
                                                      sep = ""))
# [Near from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims),
                                   mus_1b  = rep(0,n_sims), 
                                   sigmas_1b = runif(n_sims, 1, 20),
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2b  = rep(1,n_sims), 
                                   sigmas_2b = runif(n_sims, 1, 20))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     fig_name = paste(fig_basename, "_2b_esize_",
                                                      "contest_sigma_near_zero.tiff", 
                                                      sep = ""))


# Contest 3) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means and mean difference in means [Free form]
# [Far from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims, 1,1), 
                                   sigmas_1a = runif(n_sims, 1,1), 
                                   mus_1d  = runif(n_sims, 5,10), 
                                   sigmas_1d = runif(n_sims, 5,50),
                                   mus_2a  = runif(n_sims, 1,1),  
                                   sigmas_2a = runif(n_sims, 1, 1),
                                   mus_2d  = runif(n_sims, 5,10), 
                                   sigmas_2d = runif(n_sims, 5,50), ) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                     y_ax_str = "mu[d]",
                                     fig_name = paste(fig_basename, "_3a_esize_",
                                                      "contest_mu_free_near_zero.tiff", 
                                                      sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     fig_name = paste(fig_basename, "_3a_esize_",
                                                      "contest_sigma_free_near_zero.tiff", 
                                                      sep = ""))
# [Near from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,1,1), 
                                   sigmas_1a = runif(n_sims,1,1), 
                                   mus_1d  = runif(n_sims,-1,1), 
                                   sigmas_1d = runif(n_sims,1,5),
                                   mus_2a  = runif(n_sims,1,1),  
                                   sigmas_2a = runif(n_sims,1,1),
                                   mus_2d  = runif(n_sims,-1,1), 
                                   sigmas_2d = runif(n_sims,1,5) ) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_md2gtmd1", 
                                     y_ax_str = "mu[d]",
                                     fig_name = paste(fig_basename, "_3b_esize_",
                                                      "contest_mu_free_near_zero.tiff", 
                                                      sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "sigma[d]",
                                     fig_name = paste(fig_basename, "_3b_esize_",
                                                      "contest_sigma_free_near_zero.tiff", 
                                                      sep = ""))




###############################################################################
#
# Realative Error
#
##############################################################################

# Contest 4) Quantify Error rate with each metric predicting experiment with
# Lower relative mean difference in means
# [ Far from zero ]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,10,40), 
                                   sigmas_1d = rep(0,n_sims),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,50, 200), 
                                   sigmas_2d = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "r*mu[d]",
                                     fig_name = paste(fig_basename, "_4b_esize_",
                                                      "contest_rmu_far_zero.tiff",
                                                      sep = ""))
# [Near zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,-5,10), 
                                   sigmas_1d = rep(0,n_sims),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,-25, 50), 
                                   sigmas_2d = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "r*mu[d]",
                                     fig_name = paste(fig_basename, "_4b_esize_",
                                                      "contest_rmu_near_zero.tiff",
                                                      sep = ""))

# Contest 5) Quantify Error rate with each metric predicting experiment with
# Lower realtive STD of difference in means
# [ Far from zero ]
#
#------------------------------------------------------------------------------
source("R/effect_size_contests.R")
df_init <- generateExperiment_Data(n_samples=100, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(10,n_sims), 
                                    sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(10,n_sims), 
                                    sigmas_1d = runif(n_sims, 3, 5),
                                   mus_2a  = rep(40,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(10,n_sims), 
                                    sigmas_2d = runif(n_sims, 12, 20))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                      y_ax_str = "r*sigma[d]",
                      fig_name = paste(fig_basename, "_5a_esize_contest_rsigma_",
                                       "far_zero.tiff", sep = ""))

# [u_d near zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1, n_sims),
                                   mus_1d  = runif(n_sims, -1, 1), 
                                   sigmas_1d = runif(n_sims, 1, 20),
                                   mus_2a  = rep(40,n_sims), 
                                   sigmas_2a = rep(1, n_sims),
                                   mus_2d  =  runif(n_sims, -1, 1), 
                                   sigmas_2d = runif(n_sims, 4, 80))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                                     y_ax_str = "r*sigma[d]",
                                     fig_name = paste(fig_basename, "_5b_esize_",
                                                      "contest_rsigma_near_zero.tiff", 
                                                      sep = ""))



# Contest 6) Quantify Error rate with each metric predicting experiment with
# Lower rel. STD, and rel. mean, of difference in means [Free form]
# [Far from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims, 10,10), 
                                   sigmas_1a = runif(n_sims, 1, 1), 
                                   mus_1d  = runif(n_sims, 10, 25), 
                                   sigmas_1d = runif(n_sims, 2, 20),
                                   
                                   mus_2a  = runif(n_sims, 50,50),  
                                   sigmas_2a = runif(n_sims, 1, 1),
                                   mus_2d  = runif(n_sims, 50, 125), 
                                   sigmas_2d = runif(n_sims, 10, 100), ) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "r*mu[d]",
                                     fig_name = paste(fig_basename, "_6a_esize_",
                                                      "contest_rmu_free_near_zero.tiff", 
                                                      sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                                     y_ax_str = "r*sigma[d]",
                                     fig_name = paste(fig_basename, "_6a_esize_",
                                                      "contest_rsigma_free_near_zero.tiff", 
                                                      sep = ""))
# [Near from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,1,1), 
                                   sigmas_1a = runif(n_sims,1,1), 
                                   mus_1d  = runif(n_sims,-1,1), 
                                   sigmas_1d = runif(n_sims,1,5),
                                   mus_2a  = runif(n_sims,1,1),  
                                   sigmas_2a = runif(n_sims,1,1),
                                   mus_2d  = runif(n_sims,-1,1), 
                                   sigmas_2d = runif(n_sims,1,5), ) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_md2gtmd1", 
                                     y_ax_str = "r*mu[d]",
                                     fig_name = paste(fig_basename, "_6b_esize_",
                                                      "contest_rmu_free_near_zero.tiff", 
                                                      sep = ""))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma_md2gtmd1", 
                                     y_ax_str = "r*sigma[d]",
                                     fig_name = paste(fig_basename, "_6b_esize_",
                                                      "contest_rsigma_free_near_zero.tiff", 
                                                      sep = ""))











# Contest 7) which metric is best at discerning experiments below a threshold 
# difference in means
#
#------------------------------------------------------------------------------




# Contest 8) which metric is best at discerning experiments below a threshold 
# relative difference in means
#
#------------------------------------------------------------------------------




