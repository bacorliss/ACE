# sfigure_relative_agreement_contests




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
fig_num = "10" 
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


df_relative_null = list();



###############################################################################
#
# Relative Error
#
##############################################################################

# Contest 1) Lower mu_dm
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rmudm_2gt1"
fig_name = paste("F", fig_num, "_1_esize_contest_rmu_null", sep = "")
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
df_relative_null[[1]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)







# Contest 2) Lower sigma_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rsigmad_2gt1"
fig_name = paste("F", fig_num, "_2_esize_contest_rsigma_null", sep = "")
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
df_relative_null[[2]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "r*sigma[pool]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)





# Contest 3) Lower df_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
n1 <- runif(n_sims, 6, 60)
n2 <- runif(n_sims, 6, 60)
gt_colnames = "is_dfdm_2lt1"
fig_name = paste("F", fig_num, "_3_esize_contest_df_null", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed,  
                                   mus_1a  = 100, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 2, 
                                   sigmas_1ao = 5,
                                   
                                   mus_2a  = 500, 
                                   sigmas_2a = 1,
                                   mus_2ao  =  10, 
                                   sigmas_2ao = 29,
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
df_relative_null[[3]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                                   y_ax_str = "df[pool]", comp_dir ="Greater",
                                                   include_bf = include_bf, parallel_sims = parallel_sims,
                                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                                   fig_path = fig_path)




# Contest 4) Lower mu_dm, sigma_pool, df_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
gt_colnames = c("is_rmudm_2gt1","is_rsigmad_2gt1", "is_dfdm_2lt1")
fig_name = paste("F", fig_num, "_4_esize_contest_free_null", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 10, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,1,3), 
                                   sigmas_1ao = runif(n_sims,7,19),
                                   
                                   mus_2a  = 40,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,4,12),  
                                   sigmas_2ao = runif(n_sims,31,79),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
df_relative_null[[4]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                        y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                        include_bf = include_bf, parallel_sims = parallel_sims,
                                        fig_name = paste(fig_name, "_rmu.tiff",sep = ""),
                                        fig_path = fig_path)
df_relative_null[[5]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                           y_ax_str = "r*sigma[pool]",
                                           include_bf = include_bf, parallel_sims = parallel_sims,
                                           fig_name = paste(fig_name, "_rsigma.tiff",sep = ""),
                                           fig_path = fig_path)
df_relative_null[[6]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[3], 
                                                   y_ax_str = "df[pool]", comp_dir ="Greater",
                                                   include_bf = include_bf, parallel_sims = parallel_sims,
                                                   fig_name = paste(fig_name, "_df.tiff",sep = ""),
                                                   fig_path = fig_path)



save(df_relative_null, file = "temp/df_relative_null.RDS")


