

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
# p_load(broom)
p_load(tidyr)
p_load(dplyr)
# p_load(boot)
p_load(readr)
p_load(gplots)
# User defined libraries
source("R/aces.R")
source("R/agreement_contests.R")

# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
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
set.seed(rand.seed)
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_1_esize_contest_mu_near_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 20, 
                                   sigmas_1a = 10, 
                                   mus_1ao  = runif(n_sims, .5, 3), 
                                   sigmas_1ao = 10,
                                   mus_2a  = 200, 
                                   sigmas_2a = 10,
                                   mus_2ao  = runif(n_sims, 2, 4.5), 
                                   sigmas_2ao = 10,
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs,
                                   alpha_1 = 0.05, alpha_2 = 0.05,
                                   
                                   toggle_sign_rmu_d_hold_sigma = FALSE,
                                   toggle_sign_mean_ab = TRUE,
                                   switch_group_ab = TRUE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = TRUE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = FALSE,
                                   switch_n_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), 
                                   fig_path = fig_path,gt_colnames=gt_colnames)
df_unscaled_null[[1]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims, #parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)



# Contest 2) Lower sigma_d
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_sigmad_1ldt2" 
fig_name = paste("F", fig_num, "_2_esize_", "contest_sigma_near_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 10, 
                                   sigmas_1a = runif(n_sims, 1, 2),
                                   mus_1ao  = seq(.2,1.,length.out = n_sims), 
                                   sigmas_1ao = runif(n_sims, 3, 5),
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = runif(n_sims, 1, 2),
                                   mus_2ao  = seq(.2,1.,length.out = n_sims),
                                   sigmas_2ao = runif(n_sims, 4, 8),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs, 
                                   alpha_1 = 0.05, alpha_2 = 0.05,

                                   toggle_sign_rmu_d_hold_sigma = FALSE,
                                   toggle_sign_mean_ab = TRUE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_sigma_ab_12 = TRUE,
                                   switch_alpha_12 = FALSE,
                                   switch_n_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   fig_path = fig_path,gt_colnames=gt_colnames)  
df_unscaled_null[[2]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "sigma[D]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# Contest 3) Higher df_pool
# [Near from zero]
#------------------------------------------------------------------------------
n1 <- round(runif(n_sims, 10, 20))
n2 <- round(runif(n_sims, 20, 30))
set.seed(rand.seed)
gt_colnames = "is_dfdm_1ldt2"
fig_name = paste("F", fig_num, "_3_esize_", "contest_df_near_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 10, 
                                   sigmas_1a = 1,
                                   mus_1ao  = seq(0.05,.6,length.out = n_sims),
                                   sigmas_1ao = 1,
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = 1,
                                   mus_2ao  = seq(0.05,.6,length.out = n_sims),
                                   sigmas_2ao = 1,
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2, 
                                   alpha_1 = 0.05, alpha_2 = 0.05,
                                   
                                   toggle_sign_rmu_d_hold_sigma = FALSE,
                                   toggle_sign_mean_ab = TRUE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = FALSE,
                                   switch_n_12 = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), 
                                   fig_path = fig_path, gt_colnames=gt_colnames)  
df_unscaled_null[[3]] <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "df[D]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)



# Contest 4) Higher alpha
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_alpha_1ldt2"
fig_name = paste("F", fig_num, "_4_esize_", "contest_alpha_near_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 10, 
                                   sigmas_1a = .1,
                                   mus_1ao  = seq(.05,0.3,length.out = n_sims), 
                                   sigmas_1ao = 1,#seq(2,12,length.out = n_sims),
                                   
                                   mus_2a  = 100, 
                                   sigmas_2a = .1,
                                   mus_2ao  = seq(.05,0.3,length.out = n_sims),
                                   sigmas_2ao = 1,#seq(2,12,length.out = n_sims),
                                   
                                   n_1a = 30, n_1b = 30,
                                   n_2a = 30, n_2b = 30, 
                                   alpha_1 = 0.05/1,
                                   alpha_2 = 0.05/runif(n_sims, 5, 8),
                                   

                                   toggle_sign_rmu_d_hold_sigma = FALSE,
                                   toggle_sign_mean_ab = TRUE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = FALSE,
                                   switch_sigma_ab_12 = FALSE,
                                   switch_alpha_12 = TRUE,
                                   switch_n_12 = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), 
                                   fig_path = fig_path, gt_colnames=gt_colnames)  
df_unscaled_null[[4]] <-
  process_agreement_contest(df_init, gt_colname = gt_colnames, y_ax_str = "alpha[DM]",
                            include_bf = include_bf, parallel_sims = parallel_sims, 
                            fig_name = paste(fig_name, ".tiff",sep = ""),
                            fig_path = fig_path)







# Contest 5-8) Lower mu_dm, sigma_pool, df_pool
# [Near from zero]
#
#------------------------------------------------------------------------------
n1 <- round(runif(n_sims, 6, 10))
n2 <- round(runif(n_sims, 10, 14))
set.seed(rand.seed+1)
gt_colnames = c("is_mudm_1ldt2","is_sigmad_1ldt2", "is_dfdm_1ldt2","is_alpha_1ldt2")
fig_name = paste("F", fig_num, "_5_esize_contest_free_near_zero", sep = "")
df_init <- generate_population_configs(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 40, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,1, 4.5), 
                                   sigmas_1ao = runif(n_sims,10,14),
                                   
                                   mus_2a  = 1000,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,5.5, 8), 
                                   sigmas_2ao = runif(n_sims,13,16),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   alpha_1 = 0.05/runif(n_sims, 1, 2),
                                   alpha_2 = 0.05/runif(n_sims, 5, 10),
                                   
                                   toggle_sign_rmu_d_hold_sigma = FALSE,
                                   toggle_sign_mean_ab = TRUE,
                                   switch_group_ab = FALSE,
                                   switch_mu_ab_12 = FALSE,
                                   switch_mu_d_12 = TRUE,
                                   switch_sigma_ab_12 = TRUE,
                                   switch_alpha_12 = TRUE,
                                   switch_n_12 = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), 
                                   fig_path = fig_path, gt_colnames=gt_colnames)
df_unscaled_null[[5]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[1], y_ax_str = "abs(~mu[DM]*phantom(.))",
                            include_bf = include_bf, parallel_sims = parallel_sims,
                            fig_name = paste(fig_name, "_mu.tiff",sep = ""),
                            fig_path = fig_path)
df_unscaled_null[[6]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[2],  y_ax_str = "sigma[D]",
                            include_bf = include_bf, parallel_sims = parallel_sims,
                            fig_name = paste(fig_name, "_sigma.tiff",sep = ""),
                            fig_path = fig_path)
df_unscaled_null[[7]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[3], y_ax_str = "df[D]",
                            include_bf = include_bf, parallel_sims = parallel_sims,
                            fig_name = paste(fig_name, "_df.tiff",sep = ""),
                            fig_path = fig_path)
df_unscaled_null[[8]] <- 
  process_agreement_contest(df_init, gt_colname = gt_colnames[4], y_ax_str = "alpha[DM]",
                            include_bf = include_bf, parallel_sims = parallel_sims, 
                            fig_name = paste(fig_name, "_alpha.tiff",sep = ""),
                            fig_path = fig_path)




dir.create(paste(base_dir, "/temp/",sep=""),recursive = TRUE,showWarnings = FALSE)
save(df_unscaled_null, file = paste(base_dir, "/temp/df_unscaled_null.RDS",sep=""))



