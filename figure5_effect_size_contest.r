

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

source("R/effect_size_contests.R")


# Contest 1-1) Quantify Error rate with each metric discerining experiment
# with lower mean difference in means 
# [Far from zero]
#

#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas2_1a = rep(1,n_sims), 
                                   mus_1offset  = runif(n_sims,5,10), 
                                    sigmas2_1offset = rep(0,n_sims),
                                   mus_2a  = rep(1,n_sims), 
                                    sigmas2_2a = rep(1,n_sims),
                                   mus_2offset  = runif(n_sims,5,10), 
                                    sigmas2_2offset = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_d2gtd1", 
                                 y_ax_label = expression(Error~Rate~(Lower~mu[d])),
                                 out_path="temp/EffectSizeContest_mean_shift.rds", 
                                 fig_name = paste(fig_basename, "_1a_esize_contest_mu_far_zero.tiff", sep = ""))
# [Near zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas2_1a = rep(1,n_sims), 
                                   mus_1offset  = runif(n_sims,-1,1), 
                                    sigmas2_1offset = rep(0,n_sims),
                                   mus_2a  = rep(1,n_sims), 
                                    sigmas2_2a = rep(1,n_sims),
                                   mus_2offset  = runif(n_sims,-1,1), 
                                    sigmas2_2offset = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_mud_d2gtd1", 
                                     y_ax_label = expression(Error~Rate~(Lower~mu[d])),
                                     out_path="temp/EffectSizeContest_mean_shift.rds", 
                                     fig_name = paste(fig_basename, "_1b_esize_contest_mu_near_zero.tiff", sep = ""))


4
# Contest 2) which metric is better at discerining exp. with lower mean of
# [Near zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas2_1a = rep(1,n_sims), 
                                   mus_1offset  = runif(n_sims,-5,5), 
                                   sigmas2_1offset = rep(0,n_sims),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas2_2a = rep(1,n_sims),
                                   mus_2offset  = runif(n_sims,-25, 50), 
                                   sigmas2_2offset = rep(0,n_sims)) 
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rmud_d2gtd1", 
                                     y_ax_label = expression(Error~Rate~(Lower~r*mu[d])),
                                     out_path="temp/EffectSizeContest_rmean_shift.rds", 
                                     fig_name = paste(fig_basename, "_2b_esize_contest_rmu_near_zero.tiff", sep = ""))




# Contest 3) which metric is better at discerining exp. with lower 
# STD of difference in means
# [Far from zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                       mus_1a  = rep(1,n_sims), 
                                        sigmas2_1a = rep(1,n_sims),
                                       mus_1offset  = rep(99,n_sims), 
                                        sigmas2_1offset = runif(n_sims, 1, 20),
                                       mus_2a  = rep(1,n_sims), 
                                        sigmas2_2a = rep(1,n_sims),
                                       mus_2offset  = rep(99,n_sims), 
                                        sigmas2_2offset = runif(n_sims, 1, 20))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma2_d2gtd1", 
                      y_ax_label = expression(Error~Rate~(Lower~sigma[d])),
                      out_path="temp/EffectSizeContest_sigma.rds", 
                      fig_name = paste(fig_basename, "_3a_esize_contest_sigma_far_zero.tiff", sep = ""))
# [Near zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas2_1a = runif(n_sims, 1, 20),
                                   mus_1b  = rep(1,n_sims), 
                                   sigmas2_1b = runif(n_sims, 1, 20),
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas2_2a = runif(n_sims, 1, 20),
                                   mus_2b  = rep(1,n_sims), 
                                   sigmas2_2b = runif(n_sims, 1, 20))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_sigma2_d2gtd1", 
                                     y_ax_label = expression(Error~Rate~(Lower~sigma[d])),
                                     out_path="temp/EffectSizeContest_sigma.rds", 
                                     fig_name = paste(fig_basename, "_3b_esize_contest_sigma_near_zero.tiff", sep = ""))







# Contest 4) which metric is better at discerining exp. with lower 
# REL STD of difference in means
# [u_d Far zero]
#
#------------------------------------------------------------------------------
source("R/effect_size_contests.R")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas2_1a = rep(1,n_sims),
                                   mus_1offset  = rep(9,n_sims), 
                                    sigmas2_1offset = runif(n_sims, 1, 2),
                                   mus_2a  = rep(1,n_sims), 
                                    sigmas2_2a = rep(1,n_sims),
                                   mus_2offset  = rep(99,n_sims), 
                                    sigmas2_2offset = runif(n_sims, 10, 20))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigma2_d2gtd1", 
                      y_ax_label = expression(Error~Rate~(Lower~r*sigma[d])),
                      out_path="temp/EffectSizeContest_rsigma.rds", 
                      fig_name = paste(fig_basename, "_4a_esize_contest_rsigma_far_zero.tiff", sep = ""))


# [u_d near zero]
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1, n_sims), 
                                   sigmas2_1a = rep(1, n_sims), 
                                   mus_1b  = rep(10, n_sims), 
                                   sigmas2_1b = runif(n_sims, 2, 10),
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas2_2a = rep(4, n_sims),
                                   mus_2b  = rep(10, n_sims), 
                                   sigmas2_2b = runif(n_sims, 8, 40))
all_dfs <- process_esize_simulations(df_init, gt_colname = "is_rsigmad_d2gtd1", 
                                     y_ax_label = expression(Error~Rate~(Lower~r*sigma[d])),
                                     out_path="temp/EffectSizeContest_rsigma.rds", 
                                     fig_name = paste(fig_basename, "_4b_esize_contest_rsigma_near_zero.tiff", sep = ""))














# Contest 5) which metric is best at discerning experiments below a threshold 
# difference in means
#
#------------------------------------------------------------------------------




# Contest 6) which metric is best at discerning experiments below a threshold 
# relative difference in means
#
#------------------------------------------------------------------------------




