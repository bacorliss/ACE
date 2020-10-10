# sfigure_relative_agreement_contests




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


# Figure parameters
fig_num = "5" 
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
rscale_contest_path = paste("figure/SF", fig_num, "/SF", fig_num,"_rscale_contest_results.csv",sep="")


###############################################################################
#
# Relative Error
#
##############################################################################

# Contest 4) Quantify Error rate with each metric predicting experiment with
# Lower relative mean difference in means
# [Near zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rmud_md2gtmd1"
fig_name = paste("F", fig_num, "_1a_esize_contest_rmu_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,1,8), 
                                   sigmas_1ao = 24,
                                   
                                   mus_2a  = 1000, 
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,10,80), 
                                   sigmas_2ao = 240,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames)  
dfs_4a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)




# [ Far from zero ]
#
#-----------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rmud_md2gtmd1"
fig_name = paste("F", fig_num, "_1b_esize_contest_rmu_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,2,5), 
                                   sigmas_1ao = 3,
                                   
                                   mus_2a  = 500, 
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,10,25), 
                                   sigmas_2ao = 19,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)  
dfs_4b <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# Contest 5) Quantify Error rate with each metric predicting experiment with
# Lower relative STD of difference in means
# [u_d near zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rsigma_md2gtmd1"
fig_name = paste("F", fig_num, "_2a_esize_contest_rsigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed,  
                                   mus_1a  = 100, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 2, 
                                   sigmas_1ao = runif(n_sims, 5, 39),
                                   
                                   mus_2a  = 500, 
                                   sigmas_2a = 1,
                                   mus_2ao  =  10, 
                                   sigmas_2ao = runif(n_sims, 29, 199),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
dfs_5a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "r*sigma[DM]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)
# [ Far from zero ]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rsigma_md2gtmd1"
fig_name = paste("F", fig_num, "_2b_esize_contest_rsigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 20, 
                                   sigmas_1ao = runif(n_sims, 9, 49),
                                   
                                   mus_2a  = 400, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 80, 
                                   sigmas_2ao = runif(n_sims, 39, 199),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
dfs_5b <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                                    y_ax_str = "r*sigma[DM]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)



# Contest 6) Quantify Error rate with each metric predicting experiment with
# Lower rel. STD, and rel. mean, of difference in means [Free form]
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = c("is_rmud_md2gtmd1","is_rsigma_md2gtmd1")
fig_name = paste("F", fig_num, "_3a_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 10, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,1,3), 
                                   sigmas_1ao = runif(n_sims,7,19),
                                   
                                   mus_2a  = 40,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,4,12),  
                                   sigmas_2ao = runif(n_sims,31,79),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
dfs_6a_rmu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                        y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                        include_bf = include_bf, parallel_sims = parallel_sims,
                                        fig_name = paste(fig_name, "_rmu.tiff",sep = ""),
                                        fig_path = fig_path)
dfs_6a_rsigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                           y_ax_str = "r*sigma[DM]",
                                           include_bf = include_bf, parallel_sims = parallel_sims,
                                           fig_name = paste(fig_name, "_rsigma.tiff",sep = ""),
                                           fig_path = fig_path)

# [Far from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = c("is_rmud_md2gtmd1","is_rsigma_md2gtmd1")
fig_name = paste("F", fig_num, "_3b_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 20, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,4,8), 
                                   sigmas_1ao = runif(n_sims,3,11),
                                   
                                   mus_2a  = 80,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,16,32),  
                                   sigmas_2ao = runif(n_sims,15,47),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
dfs_6b_rmu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                        y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                        include_bf = include_bf, parallel_sims = parallel_sims,
                                        fig_name = paste(fig_name, "_rmu.tiff",sep = ""),
                                        fig_path = fig_path)
dfs_6b_rsigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                           y_ax_str = "r*sigma[DM]",
                                           include_bf = include_bf, parallel_sims = parallel_sims,
                                           fig_name = paste(fig_name, "_rsigma.tiff",sep = ""),
                                           fig_path = fig_path)




# Save and load effect size contest data
save.image(file=paste("temp/F", fig_num, "_effect_size_contest_results.rds"))
load(paste("temp/F", fig_num, "_effect_size_contest_results.rds"))

