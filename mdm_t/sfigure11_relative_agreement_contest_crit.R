


#' Relative agreement contests, critical region
#' Calculates comparison error for candidate statistics in determining which of 
#' two results have higher relative agreement. Each investigation varies single agreement 
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
source("R/aces.R")
source("R/agreement_contests.R")
base_dir = "mdm_t"


# Figure parameters
#-------------------------------------------------------------------------------
fig_num = "11" 
fig_path = paste(base_dir, "/figure/SF",fig_num, "/",sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Simulation parameters
# 
#-------------------------------------------------------------------------------
n_sims = 1e3
n_samples = 1e2
n_obs = 50
rand.seed = 1

parallel_sims = TRUE
include_bf = TRUE
rscale_contest_path = paste(base_dir, "/figure/SF", fig_num, "/SF", fig_num,
                            "_rscale_contest_results.csv",sep="")
df_relative_crit = list();



# Contest 1: lower rmu_dm
#
#-----------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rmudm_1ldt2"
fig_name = paste("F", fig_num, "_1_esize_contest_rmu_crit", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1, 
                                   rmus_1d  = runif(n_sims, 0.01, 0.018), 
                                   rsigmas_1d = 0.035,
                                   
                                   mus_2a  = 1000, 
                                   sigmas_2a = 1,
                                   rmus_2d  = runif(n_sims, 0.017, 0.025), 
                                   rsigmas_2d = 0.035,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   alpha_1 = 0.05, alpha_2 = 0.05,
                                   
                                   toggle_sign_rmu_d_hold_sigma = TRUE,
                                   toggle_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_rmu_d_12_hold_rsigma = TRUE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = FALSE,
                                   switch_n_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)  
df_relative_crit[[1]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)




# Contest 2: lower rsigma_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rsigmad_1ldt2"
fig_name = paste("F", fig_num, "_2_esize_contest_rsigma_crit", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1,
                                   rmus_1d  = 0.15, 
                                   rsigmas_1d = runif(n_sims, 0.25, 0.35),
                                   
                                   mus_2a  = 500, 
                                   sigmas_2a = 1,
                                   rmus_2d  = 0.15, 
                                   rsigmas_2d = runif(n_sims, 0.35, 0.5),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   alpha_1 = 0.05, alpha_2 = 0.05,
                                   
                                   toggle_sign_rmu_d_hold_sigma = TRUE,
                                   toggle_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_rmu_d_12_hold_rsigma = FALSE,
                                   switch_rsigma_ab_12_hold_sigma_a = TRUE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = FALSE,
                                   switch_n_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
df_relative_crit[[2]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "r*sigma[D]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# Contest 3: lower df_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
n1 <- runif(n_sims, 6, 15)
n2 <- runif(n_sims, 15, 30)
gt_colnames = "is_dfdm_1ldt2"
fig_name = paste("F", fig_num, "_3_esize_contest_df_crit", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1,
                                   rmus_1d  = 0.45, 
                                   rsigmas_1d = runif(n_sims, 0.25, 0.45),
                                   
                                   mus_2a  = 400, 
                                   sigmas_2a = 1,
                                   rmus_2d  = 0.45, 
                                   rsigmas_2d = runif(n_sims, 0.25, 0.45),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   toggle_sign_rmu_d_hold_sigma = TRUE,
                                   toggle_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_rmu_d_12_hold_rsigma = FALSE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = FALSE,
                                   switch_n_12 = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
df_relative_crit[[3]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                                   y_ax_str = "df[D]", 
                                                   include_bf = include_bf, parallel_sims = parallel_sims,
                                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                                   fig_path = fig_path)



# Contest 4) Lower alpha_dm
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_alpha_1ldt2"
fig_name = paste("F", fig_num, "_4_esize_contest_alpha_crit", sep = "")
mus_1ao = seq(0.08,0.01, length.out = n_sims)
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1, 
                                   rmus_1d  = seq(0.02, 0.06, length.out = n_sims),
                                   rsigmas_1d = 0.05,
                                   
                                   mus_2a  = 1000, 
                                   sigmas_2a = 1,
                                   rmus_2d  = seq(0.02, 0.06, length.out = n_sims),
                                   rsigmas_2d = 0.05,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   alpha_1 = 0.05/runif(n_sims, 1, 2),
                                   alpha_2 = 0.05/runif(n_sims, 5, 10),
                                   
                                   toggle_sign_rmu_d_hold_sigma = TRUE,
                                   toggle_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_rmu_d_12_hold_rsigma = FALSE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = TRUE,
                                   switch_n_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames)  
df_relative_crit[[4]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                                   y_ax_str = "alpha[DM]",
                                                   include_bf = include_bf, parallel_sims = parallel_sims,
                                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                                   fig_path = fig_path)



# Contest 6) Lower mu_dm, sigma_pool, df_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed+5)
n1 <- runif(n_sims, 6, 25)
n2 <- runif(n_sims, 30, 50)
gt_colnames = c("is_rmudm_1ldt2","is_rsigmad_1ldt2", "is_dfdm_1ldt2","is_alpha_1ldt2")
fig_name = paste("F", fig_num, "_5_esize_contest_free_crit", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 10, 
                                   sigmas_1a = 1, 
                                   rmus_1d  = runif(n_sims, 0.36, 0.45), 
                                   rsigmas_1d = runif(n_sims, 0.25, 0.3),
                                   
                                   mus_2a  = 50,  
                                   sigmas_2a = 1,
                                   rmus_2d  = runif(n_sims, 0.45, 0.5), 
                                   rsigmas_2d = runif(n_sims, 0.4, 0.46),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   alpha_1 = 0.05/runif(n_sims, 1, 2),
                                   alpha_2 = 0.05/runif(n_sims, 5, 10),
                                   
                                   toggle_sign_rmu_d_hold_sigma = TRUE,
                                   toggle_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_rmu_d_12_hold_rsigma = TRUE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_rsigma_ab_12_hold_sigma_a = TRUE,
                                   switch_alpha_12 = TRUE,
                                   switch_n_12 = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
df_relative_crit[[5]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[1], y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                            include_bf = include_bf, parallel_sims = parallel_sims,
                            fig_name = paste(fig_name, "_rmu.tiff",sep = ""),
                            fig_path = fig_path)
df_relative_crit[[6]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[2], y_ax_str = "r*sigma[D]",
                            include_bf = include_bf, parallel_sims = parallel_sims,
                            fig_name = paste(fig_name, "_rsigma.tiff",sep = ""),
                            fig_path = fig_path)
df_relative_crit[[7]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[3], y_ax_str = "df[D]", 
                            include_bf = include_bf, parallel_sims = parallel_sims,
                            fig_name = paste(fig_name, "_df.tiff",sep = ""),
                            fig_path = fig_path)
df_relative_crit[[8]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[4], y_ax_str = "alpha[DM]",
                            include_bf = include_bf, parallel_sims = parallel_sims, 
                            fig_name = paste(fig_name, "_alpha.tiff",sep = ""),
                            fig_path = fig_path)










# Output results
dir.create(paste(base_dir, "/temp/",sep=""),recursive = TRUE,showWarnings = FALSE)
save(df_relative_crit, file = paste(base_dir, "/temp/df_relative_crit.RDS",sep=""))