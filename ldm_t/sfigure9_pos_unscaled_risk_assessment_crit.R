


#' Unscaled agreement contests, critical region
#' Calculates comparison error for candidate statistics in determining which of 
#' two results have higher unscaled agreement. Each investigation varies single agreement 
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
source("R/strength_risk_assessment.R")

# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "ldm_t"
fig_num = "9" 
fig_path = paste(base_dir, "/figure/SF",fig_num, "/",sep="")
dir.create(file.path(getwd(), fig_path), showWarnings = FALSE, recursive = TRUE)

# Simulation parameters
#-------------------------------------------------------------------------------
n_sims = 1e3
if (!exists("n_samples",mode = "numeric")) {n_samples = 1e3}
if (!exists("use_pseudo_samples")) {use_pseudo_samples = TRUE}
n_obs = 50
rand.seed = 1
parallel_sims = TRUE
scale_contest_path = paste(base_dir, "/figure/SF", fig_num, "/SF", fig_num,
                           "_scale_contest_results.csv",sep="")

df_unscaled_pos = list();

# Contest 1) Lower mu_d
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_mudm_1hest2"
fig_name = paste("F", fig_num, "_1_esize_contest_mu_far_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                       mus_1a  = 10, 
                                       sigmas_1a = 2, 
                                       mus_1ao  = runif(n_sims, 1, 3), 
                                       sigmas_1ao = 2,
                                       mus_2a  = 300, 
                                       sigmas_2a = 2,
                                       mus_2ao  = runif(n_sims, 3, 5), 
                                       sigmas_2ao = 2,
                                       
                                       n_1a = n_obs, n_1b = n_obs,
                                       n_2a = n_obs, n_2b = n_obs,
                                       alpha_1 = 0.05, alpha_2 = 0.05,
                                       
                                       toggle_sign_rmu_d_hold_rsigma = FALSE,
                                       toggle_sign_mean_ab = FALSE,
                                       switch_group_ab = FALSE,
                                       switch_mu_ab_12 = FALSE,
                                       switch_mu_d_12 = TRUE,
                                       switch_sigma_ab_12 = FALSE,
                                       switch_alpha_12 = FALSE,
                                       switch_n_12 = FALSE,
                                       fig_name = paste(fig_name, ".tiff",sep = ""), 
                                       fig_path = fig_path,gt_colnames=gt_colnames,
                                       strength = "hest" )
df_unscaled_pos[[1]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames, 
                           measure_pretty_str = "abs(~mu[DM]*phantom(.))",
                           parallel_sims = parallel_sims, #parallel_sims,
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)





# Contest 2) Lower sigma_pooled
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_sigmad_1hest2" 
fig_name = paste("F", fig_num, "_2_esize_", "_2_esize_contest_sigma_far_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                       mus_1a  = 10, 
                                       sigmas_1a = 0.5,
                                       mus_1ao  = seq(0.7, 1, n_sims), 
                                       sigmas_1ao = runif(n_sims, 0.4, 1),
                                       
                                       mus_2a  = 100, 
                                       sigmas_2a = 0.5,
                                       mus_2ao  = seq(0.7, 1, n_sims),
                                       sigmas_2ao = runif(n_sims, 0.9, 1.8),
                                       
                                       n_1a = n_obs, n_1b = n_obs,
                                       n_2a = n_obs, n_2b = n_obs, 
                                       alpha_1 = 0.05, alpha_2 = 0.05,
                                       
                                       toggle_sign_rmu_d_hold_rsigma = FALSE,
                                       toggle_sign_mean_ab = FALSE,
                                       switch_group_ab = FALSE,
                                       switch_mu_ab_12 = FALSE,
                                       switch_mu_d_12 = FALSE,
                                       switch_sigma_ab_12 = TRUE,
                                       switch_alpha_12 = FALSE,
                                       switch_n_12 = FALSE,
                                       fig_name = paste(fig_name, ".tiff",sep = ""),
                                       fig_path = fig_path,gt_colnames=gt_colnames,
                                       strength = "hest" )  
df_unscaled_pos[[2]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames, 
                           measure_pretty_str = "sigma[D]",
                           parallel_sims = parallel_sims,
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)




# Contest 3) Higher df_pool
# [Near from zero]
source("R/aces.R")
#------------------------------------------------------------------------------
n1 <- round(runif(n_sims, 6, 12))
n2 <- round(runif(n_sims, 15, 30))
set.seed(rand.seed)
gt_colnames = "is_dfdm_1hest2"
fig_name = paste("F", fig_num, "_3_esize_contest_df_far_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                       mus_1a  = 10, 
                                       sigmas_1a = 1,
                                       mus_1ao  = seq(1.2,3,length.out = n_sims),
                                       sigmas_1ao = 1,
                                       
                                       mus_2a  = 100, 
                                       sigmas_2a = 1,
                                       mus_2ao  = seq(1.2,3,length.out = n_sims),
                                       sigmas_2ao = 1,
                                       
                                       n_1a = n1, n_1b = n1,
                                       n_2a = n2, n_2b = n2, 
                                       alpha_1 = 0.05, alpha_2 = 0.05,
                                       
                                       toggle_sign_rmu_d_hold_rsigma = FALSE,
                                       toggle_sign_mean_ab = FALSE,
                                       switch_group_ab = FALSE,
                                       switch_mu_ab_12 = FALSE,
                                       switch_mu_d_12 = FALSE,
                                       switch_sigma_ab_12 = FALSE,
                                       switch_alpha_12 = FALSE,
                                       switch_n_12 = TRUE,
                                       fig_name = paste(fig_name, ".tiff",sep = ""), 
                                       fig_path = fig_path, gt_colnames=gt_colnames,
                                       strength = "hest" ) 
df_unscaled_pos[[3]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames, 
                           measure_pretty_str = "df[D]",
                           parallel_sims = TRUE,
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)











# Contest 4) Higher alpha
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_alpha_1hest2"
fig_name = paste("F", fig_num, "_4_esize_", "contest_alpha_far_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                       mus_1a  = 10, 
                                       sigmas_1a = .1,
                                       mus_1ao  = seq(.6, 1.5,length.out = n_sims), 
                                       sigmas_1ao = 1,#seq(2,12,length.out = n_sims),
                                       
                                       mus_2a  = 500, 
                                       sigmas_2a = .1,
                                       mus_2ao  = seq(.6, 1.5,length.out = n_sims),
                                       sigmas_2ao = 1,#seq(2,12,length.out = n_sims),
                                       
                                       n_1a = 30, n_1b = 30,
                                       n_2a = 30, n_2b = 30, 
                                       alpha_1 = 0.05/runif(n_sims, 1, 2),
                                       alpha_2 = 0.05/runif(n_sims, 5, 10),
                                       
                                       
                                       toggle_sign_rmu_d_hold_rsigma = FALSE,
                                       toggle_sign_mean_ab = FALSE,
                                       switch_group_ab = FALSE,
                                       switch_mu_ab_12 = FALSE,
                                       switch_mu_d_12 = FALSE,
                                       switch_sigma_ab_12 = FALSE,
                                       switch_alpha_12 = TRUE,
                                       switch_n_12 = FALSE,
                                       fig_name = paste(fig_name, ".tiff",sep = ""), 
                                       fig_path = fig_path, gt_colnames=gt_colnames,
                                       strength = "hest" )  
df_unscaled_pos[[4]] <-
  process_strength_contest(df_init, gt_colname = gt_colnames, measure_pretty_str = "alpha[DM]",
                           parallel_sims = parallel_sims, 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)










# Contest 7) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means and mean difference in means [Free form]
# [Near from zero]
# [Far from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
n1 <- round(runif(n_sims, 5, 8))
n2 <- round(runif(n_sims, 15, 30))
gt_colnames = c("is_mudm_1hest2","is_sigmad_1hest2", "is_dfdm_1hest2","is_alpha_1hest2")
fig_name = paste("F", fig_num, "_4_esize_contest_free_far_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                       mus_1a  = 10, 
                                       sigmas_1a = 1, 
                                       mus_1ao  = runif(n_sims, 3.8, 5), 
                                       sigmas_1ao = runif(n_sims, 1, 1.75),
                                       
                                       mus_2a  = 100,
                                       sigmas_2a = 1,
                                       mus_2ao  = runif(n_sims, 4.5, 6), 
                                       sigmas_2ao = runif(n_sims, 2.25, 3),
                                       
                                       n_1a = n1, n_1b = n1,
                                       n_2a = n2, n_2b = n2,
                                       alpha_1 = 0.05/runif(n_sims, 1, 2),
                                       alpha_2 = 0.05/runif(n_sims, 5, 10),
                                       
                                       toggle_sign_rmu_d_hold_rsigma = FALSE,
                                       toggle_sign_mean_ab = FALSE,
                                       switch_group_ab = FALSE,
                                       switch_mu_ab_12 = FALSE,
                                       switch_mu_d_12 = TRUE,
                                       switch_sigma_ab_12 = TRUE,
                                       switch_alpha_12 = TRUE,
                                       switch_n_12 = TRUE,
                                       fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                       gt_colnames=gt_colnames,
                                       strength = "hest")
df_unscaled_pos[[5]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames[1], measure_pretty_str = "abs(~mu[DM]*phantom(.))",
                           parallel_sims = parallel_sims,
                           fig_name = paste(fig_name, "_mu.tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)
df_unscaled_pos[[6]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames[2], measure_pretty_str = "sigma[D]",
                           parallel_sims = parallel_sims,
                           fig_name = paste(fig_name, "_sigma.tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)
df_unscaled_pos[[7]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames[3], measure_pretty_str = "df[D]",
                           parallel_sims = parallel_sims,
                           fig_name = paste(fig_name, "_df.tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)
df_unscaled_pos[[8]] <- 
  process_strength_contest(df_init, gt_colname = gt_colnames[4], measure_pretty_str = "alpha[DM]",
                           parallel_sims = parallel_sims, 
                           fig_name = paste(fig_name, "_alpha.tiff",sep = ""),
                           fig_path = fig_path, stat_exclude_list = NULL,
                           strength = "hest", delta = 1, is_delta_relative = FALSE,
                           use_pseudo_samples = use_pseudo_samples)


dir.create(paste(base_dir, "/temp/",sep=""),recursive = TRUE,showWarnings = FALSE)
save(df_unscaled_pos, file = paste(base_dir, "/temp/df_unscaled_pos.RDS",sep=""))