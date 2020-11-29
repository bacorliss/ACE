





# Unscaled agreement null region

# Set workspace to location of this source file, since there are file dependencies
# on a relative path basis
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
base_dir = "mmd_z"

# Agreement COntests for Unscaled Error


# Figure parameters
fig_num = "9" 
fig_path = paste(base_dir, "/figure/SF",fig_num, "/",sep="")
dir.create(file.path(getwd(), fig_path), showWarnings = FALSE)
# fig_path = paste("figure/SF",fig_num, "/",sep="")
# Simulation parameters: A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_sims = 1e3
n_samples = 1e2
n_obs = 50
rand.seed = 1

parallel_sims = TRUE
include_bf = TRUE
scale_contest_path = paste(base_dir, "/figure/SF", fig_num, "/SF", fig_num,"_scale_contest_results.csv",sep="")



df_unscaled_crit = list();



# Contest 1) Lower mu_dm
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_mudm_2gt1"
fig_name = paste("F", fig_num, "_1_esize_contest_mu_far_zero", sep = "")
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
df_unscaled_crit[[1]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# COntest 2) Lower sigmad
# [Far from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_sigmad_2gt1"
fig_name = paste("F", fig_num, "_2_esize_contest_sigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 10, 
                                   sigmas_1ao = runif(n_sims, 5, 20),
                                   
                                   mus_2a  = 50, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 10,
                                   sigmas_2ao = runif(n_sims, 5, 20),
                                   
                                   n_1a = n_obs, n_1b = n_obs,
                                   n_2a = n_obs, n_2b = n_obs, 

                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
df_unscaled_crit[[2]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "sigma[pool]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)



# Contest 7) Lower df_pool
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
set.seed(rand.seed)
gt_colnames = "is_dfdm_2lt1"
fig_name = paste("F", fig_num, "_3_esize_contest_df_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 10, 
                                   sigmas_1ao = 10,
                                   
                                   mus_2a  = 50, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 10,
                                   sigmas_2ao = 10,
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2, 
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
df_unscaled_crit[[3]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                            y_ax_str = "df[pool]", comp_dir ="Greater",
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
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
gt_colnames = c("is_mudm_2gt1","is_sigmad_2gt1", "is_dfdm_2lt1")
fig_name = paste("F", fig_num, "_4_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 1, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims, 3, 6), 
                                   sigmas_1ao = runif(n_sims, 2.5, 7.5),
                                   
                                   mus_2a  = 100,
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims, 3, 6), 
                                   sigmas_2ao = runif(n_sims, 2.5, 7.5),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
df_unscaled_crit[[4]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                       y_ax_str = "abs(~mu[DM]*phantom(.))",
                                       include_bf = include_bf, parallel_sims = parallel_sims,
                                       fig_name = paste(fig_name, "_mu.tiff",sep = ""),
                                       fig_path = fig_path)
df_unscaled_crit[[5]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                          y_ax_str = "sigma[pool]",
                                          include_bf = include_bf, parallel_sims = parallel_sims,
                                          fig_name = paste(fig_name, "_sigma.tiff",sep = ""),
                                          fig_path = fig_path)
df_unscaled_crit[[6]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[3], 
                                                  y_ax_str = "df[pool]", comp_dir ="Greater",
                                                  include_bf = include_bf, parallel_sims = parallel_sims,
                                                  fig_name = paste(fig_name, "_df.tiff",sep = ""),
                                                  fig_path = fig_path)


dir.create(paste(base_dir, "/temp/",sep=""),recursive = TRUE,showWarnings = FALSE)
save(df_unscaled_crit, file = paste(base_dir, "/temp/df_unscaled_crit.RDS"))