

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



# Contest 1-1) Quantify Error rate with each metric discerining experiment
# with lower mean difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_1b  = runif(n_sims, 1, 2), 
                                   sigmas_1b = rep(1,n_sims), 
                                   mus_2b  = runif(n_sims, 1, 2), 
                                   sigmas_2b = rep(1,n_sims)) 
all_dfs <- quantifyTidyPlot_Data(df_init, gt_colname = "is_mud_d2gtd1", 
                      y_label_str = expression(Error~Rate~(Lower~mu[d])),
                      out_path="temp/EffectSizeContest_mean_shift.rds", 
                      fig_name = paste(fig_basename, "_1a effect size contest_mu.tiff", sep = ""))
  
  

# Contest 2) which metric is better at discerining exp. with lower mean of
# relative difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(50, n_sims),
                                   sigmas_1a = rep(1,n_sims),
                                   mus_2a  = rep(50, n_sims),
                                   sigmas_2a = rep(1,n_sims),
                                   mus_1b  = runif(n_sims, 10, 100),
                                   sigmas_1b = rep(1,n_sims),
                                   mus_2b  = runif(n_sims, 10, 100),
                                   sigmas_2b = rep(1,n_sims))
all_dfs <- quantifyTidyPlot_Data(df_init, gt_colname = "is_rmud_d2gtd1", 
                      y_label_str = expression(Error~Rate~(Lower~r*mu[d])),
                      out_path="temp/EffectSizeContest_rmean_shift.rds", 
                      fig_name = paste(fig_basename, "_2c effect size contest_rmu.tiff", sep = ""))


# Contest 3) which metric is better at discerining exp. with lower 
# STD of difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                       mus_1a  = rep(1,n_sims), 
                                       sigmas_1a = rep(1,n_sims), 
                                       mus_2a  = rep(5,n_sims), 
                                       sigmas_2a = rep(1,n_sims),
                                       mus_1b  = rep(1,n_sims), 
                                       sigmas_1b = runif(n_sims, 0.1, 1),
                                       mus_2b  = rep(5,n_sims), 
                                       sigmas_2b = runif(n_sims, 0.1, 1))
all_dfs <- quantifyTidyPlot_Data(df_init, gt_colname = "is_sigmad_d2gtd1", 
                      y_label_str = expression(Error~Rate~(Lower~sigma[d])),
                      out_path="temp/EffectSizeContest_sigma.rds", 
                      fig_name = paste(fig_basename, "_4d effect size contest_sigma.tiff", sep = ""))


# Contest 4) which metric is better at discerining exp. with lower 
# REL STD of difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                       mus_1a  = rep(1,n_sims), 
                                       sigmas_1a = rep(1,n_sims), 
                                       mus_2a  = rep(1,n_sims), 
                                       sigmas_2a = rep(1,n_sims),
                                       mus_1b  = rep(1,n_sims), 
                                       sigmas_1b = runif(n_sims, 1, 2),
                                       mus_2b  = rep(1,n_sims), 
                                       sigmas_2b = runif(n_sims, 9, 10))
all_dfs <- quantifyTidyPlot_Data(df_init, gt_colname = "is_rsigmad_d2gtd1", 
                      y_label_str = expression(Error~Rate~(Lower~r*sigma[d])),
                      out_path="temp/EffectSizeContest_rsigma.rds", 
                      fig_name = paste(fig_basename, "_5e effect size contest_rsigma.tiff", sep = ""))







# Contest 5) which metric is best at discerning experiments below a threshold 
# difference in means
#
#------------------------------------------------------------------------------




# Contest 6) which metric is best at discerning experiments below a threshold 
# relative difference in means
#
#------------------------------------------------------------------------------




