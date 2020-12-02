

#' Unscaled agreement contests, null region
#' Calculates comparison error for candidate statistics in determining which of 
#' two results have higherunscaled agreement. Each investigation varies single agreement 
#' parameter as independent variable (or varies multiple simultaneously). 
#' Selected independent variable(s) are used as ground truth to determine error rates.

# Load required packages
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
base_dir = "mmd_z"
fig_num = "8" 
dir.create(file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep="")), 
           showWarnings = FALSE, recursive = TRUE)
fig_path = paste(base_dir, "/figure/SF",fig_num, "/",sep="")


# Simulation parameters
#-------------------------------------------------------------------------------
n_sims = 1e3
n_samples = 1e2
n_obs = 50
rand.seed = 1
parallel_sims = TRUE
include_bf = TRUE
scale_contest_path = paste(base_dir, "/figure/SF", fig_num, "/SF", fig_num,"_scale_contest_results.csv",sep="")
df_unscaled_null = list();


# Contest 1) Lower mu_d
# [Near from zero]
#
#------------------------------------------------------------------------------
source("R/agreement_contests.R")
source("R/row_effect_sizes.R")
set.seed(rand.seed)
gt_colnames = "is_mudm_1hat2"
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
df_unscaled_null[[1]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims, #parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)



# COntest 2) Lower sigma_pooled
# [Near from zero]
#
#------------------------------------------------------------------------------

set.seed(rand.seed)
gt_colnames = "is_sigmad_1hat2" 
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
df_unscaled_null[[2]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "sigma[pool]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# COntest 3) Higher df_pool
# [Near from zero]
#
source("R/agreement_contests.R")
#------------------------------------------------------------------------------
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
set.seed(rand.seed)
gt_colnames = "is_dfdm_1hat2"
fig_name = paste("F", fig_num, "_3_esize_", "contest_df_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 0.5, 
                                   sigmas_1ao = 1,
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 0.5,
                                   sigmas_2ao = 1,
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2, 
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = FALSE,
                                   switch_exp_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)  
df_unscaled_null[[3]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "df[pool]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)




# Contest 4) Lower mu_dm, sigma_pool, df_pool
# [Near from zero]
#
#------------------------------------------------------------------------------
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
set.seed(rand.seed)
gt_colnames = c("is_mudm_1hat2","is_sigmad_1hat2", "is_dfdm_1hat2")
fig_name = paste("F", fig_num, "_4_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 2, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,0.5,10), 
                                   sigmas_1ao = runif(n_sims,24,44),
                                   
                                   mus_2a  = 50,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,0.5,10), 
                                   sigmas_2ao = runif(n_sims,24,44),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
df_unscaled_null[[4]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                       y_ax_str = "abs(~mu[DM]*phantom(.))",
                                       include_bf = include_bf, parallel_sims = parallel_sims,
                                       fig_name = paste(fig_name, "_mu.tiff",sep = ""),
                                       fig_path = fig_path)
df_unscaled_null[[5]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                          y_ax_str = "sigma[pool]",
                                          include_bf = include_bf, parallel_sims = parallel_sims,
                                          fig_name = paste(fig_name, "_sigma.tiff",sep = ""),
                                          fig_path = fig_path)
df_unscaled_null[[6]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[3], 
                                          y_ax_str = "df[pool]", 
                                          include_bf = include_bf, parallel_sims = parallel_sims,
                                          fig_name = paste(fig_name, "_df.tiff",sep = ""),
                                          fig_path = fig_path)



dir.create(paste(base_dir, "/temp/",sep=""),recursive = TRUE,showWarnings = FALSE)
save(df_unscaled_null, file = paste(base_dir, "/temp/df_unscaled_null.RDS",sep=""))



