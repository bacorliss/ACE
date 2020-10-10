

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

# Set workspace to location of this source file, since there are file dependencies
# on a relative path basis
# From #R Studio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd(getSrcDirectory()[1])
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
source("R/agreement_contests.R")

# Agreement COntests for Unscaled Error


# Figure parameters
fig_num = "4" 
dir.create(file.path(getwd(), paste("figure/SF",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/SF",fig_num, "/",sep="")
# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_sims = 1e3
n_samples = 1e2
n_obs = 50
rand.seed = 1

parallel_sims = TRUE
include_bf = TRUE
scale_contest_path = paste("figure/SF", fig_num, "/SF", fig_num,"_scale_contest_results.csv",sep="")



# Contest 1) Quantify Error rate with each metric predicting experiment with
# Lower mean difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1_esize_contest_mu_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims, 1, 8), 
                                   sigmas_1ao = 19,
                                   mus_2a  = 100, 
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims, 1, 8), 
                                   sigmas_2ao = 19,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
dfs_1a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)



# COntest 2) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_sigma_md2gtmd1"
fig_name = paste("F", fig_num, "_2_esize_", "contest_sigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 5, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 1.5, 
                                   sigmas_1ao = runif(n_sims, 3, 25),
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 1.5,
                                   sigmas_2ao = runif(n_sims, 3, 25),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs, 
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = FALSE,
                                   switch_exp_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)  
dfs_2a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "sigma[DM]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# COntest 3) Quantify Error rate with each metric predicting experiment with
# Lower degrees of freedom in the difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_sigma_md2gtmd1"
fig_name = paste("F", fig_num, "_2_esize_", "contest_sigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 1, 
                                   sigmas_1ao = 1,
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 1,
                                   sigmas_2ao = 1,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs, 
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = FALSE,
                                   switch_exp_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)  
dfs_2a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "sigma[DM]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)




# Contest 3) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means and mean difference in means [Free form]
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = c("is_mud_md2gtmd1","is_sigma_md2gtmd1")
fig_name = paste("F", fig_num, "_3_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 2, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,0.5,10), 
                                   sigmas_1ao = runif(n_sims,24,44),
                                   
                                   mus_2a  = 50,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,0.5,10), 
                                   sigmas_2ao = runif(n_sims,24,44),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
dfs_3a_mu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                       y_ax_str = "abs(~mu[DM]*phantom(.))",
                                       include_bf = include_bf, parallel_sims = parallel_sims,
                                       fig_name = paste(fig_name, "_mu.tiff",sep = ""),
                                       fig_path = fig_path)
dfs_3a_sigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                          y_ax_str = "sigma[DM]",
                                          include_bf = include_bf, parallel_sims = parallel_sims,
                                          fig_name = paste(fig_name, "_sigma.tiff",sep = ""),
                                          fig_path = fig_path)













# Contest 5) Quantify Error rate with each metric predicting experiment with
# Lower mean difference in means
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_5_esize_contest_mu_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,1,3), 
                                   sigmas_1ao = 1,
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,1,3), 
                                   sigmas_2ao = 1,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)  
dfs_1b <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# COntest 6) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means
# [Far from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_sigma_md2gtmd1"
fig_name = paste("F", fig_num, "_6_esize_contest_sigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 10, 
                                   sigmas_1ao = runif(n_sims, 5, 20),
                                   
                                   mus_2a  = 50, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 10,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs, 
                                   
                                   sigmas_2ao = runif(n_sims, 5, 20),
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
dfs_2b <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "sigma[DM]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# Contest 7) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means and mean difference in means [Free form]
# [Near from zero]
# [Far from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = c("is_mud_md2gtmd1","is_sigma_md2gtmd1")
fig_name = paste("F", fig_num, "_7_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims, 3, 6), 
                                   sigmas_1ao = runif(n_sims, 2.5, 7.5),
                                   
                                   mus_2a  = 100,
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims, 3, 6), 
                                   sigmas_2ao = runif(n_sims, 2.5, 7.5),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
dfs_3b_mu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                       y_ax_str = "abs(~mu[DM]*phantom(.))",
                                       include_bf = include_bf, parallel_sims = parallel_sims,
                                       fig_name = paste(fig_name, "_mu.tiff",sep = ""),
                                       fig_path = fig_path)
dfs_3b_sigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                          y_ax_str = "sigma[DM]",
                                          include_bf = include_bf, parallel_sims = parallel_sims,
                                          fig_name = paste(fig_name, "_sigma.tiff",sep = ""),
                                          fig_path = fig_path)