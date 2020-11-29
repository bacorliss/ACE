
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
base_dir = "mmd_z"


# Figure parameters
fig_num = "11" 
fig_path = paste(base_dir, "/figure/SF",fig_num, "/",sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

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
rscale_contest_path = paste(base_dir, "/figure/SF", fig_num, "/SF", fig_num,"_rscale_contest_results.csv",sep="")
df_relative_crit = list();


# Contest 1: lower rmu_dm
#
#-----------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rmudm_2gt1"
fig_name = paste("F", fig_num, "_1_esize_contest_rmu_crit", sep = "")
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
df_relative_crit[[1]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)




# Contest 2: lower rsigma_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_rsigmad_2gt1"
fig_name = paste("F", fig_num, "_2_esize_contest_rsigma_crit", sep = "")
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
df_relative_crit[[2]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "r*sigma[pool]",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path)


# Contest 3: lower df_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
gt_colnames = "is_dfdm_2lt1"
fig_name = paste("F", fig_num, "_3_esize_contest_df_crit", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 100, 
                                   sigmas_1a = 1,
                                   mus_1ao  = 20, 
                                   sigmas_1ao = runif(n_sims, 9, 49),
                                   
                                   mus_2a  = 400, 
                                   sigmas_2a = 1,
                                   mus_2ao  = 80, 
                                   sigmas_2ao = runif(n_sims, 39, 199),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames) 
df_relative_crit[[3]] <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                                   y_ax_str = "df[pool]", comp_dir ="Greater",
                                                   include_bf = include_bf, parallel_sims = parallel_sims,
                                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                                   fig_path = fig_path)


# Contest 4: lower rmu_dm, rsigma_pool, df_pool
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
n1 <- runif(n_sims, 6, 75)
n2 <- runif(n_sims, 6, 75)
gt_colnames = c("is_rmudm_2gt1","is_rsigmad_2gt1","is_dfdm_2lt1")
fig_name = paste("F", fig_num, "_4_esize_contest_free_crit", sep = "")
df_init <- generateExperiment_Data(n_samples=n_samples, n_sims=n_sims, rand.seed=rand.seed, 
                                   mus_1a  = 20, 
                                   sigmas_1a = 1, 
                                   mus_1ao  = runif(n_sims,4,8), 
                                   sigmas_1ao = runif(n_sims,3,11),
                                   
                                   mus_2a  = 80,  
                                   sigmas_2a = 1,
                                   mus_2ao  = runif(n_sims,16,32),  
                                   sigmas_2ao = runif(n_sims,15,47),
                                   
                                   n_1a = n1, n_1b = n1,
                                   n_2a = n2, n_2b = n2,
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""), fig_path = fig_path,
                                   gt_colnames=gt_colnames)
df_relative_crit[[4]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                        y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                        include_bf = include_bf, parallel_sims = parallel_sims,
                                        fig_name = paste(fig_name, "_rmu.tiff",sep = ""),
                                        fig_path = fig_path)
df_relative_crit[[5]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                           y_ax_str = "r*sigma[pool]",
                                           include_bf = include_bf, parallel_sims = parallel_sims,
                                           fig_name = paste(fig_name, "_rsigma.tiff",sep = ""),
                                           fig_path = fig_path)
df_relative_crit[[6]] <- process_esize_simulations(df_init, gt_colname = gt_colnames[3], 
                                           y_ax_str = "df[pool]", comp_dir ="Greater",
                                           include_bf = include_bf, parallel_sims = parallel_sims,
                                           fig_name = paste(fig_name, "_df.tiff",sep = ""),
                                           fig_path = fig_path)



# Output results
dir.create(paste(base_dir, "/temp/",sep=""),recursive = TRUE,showWarnings = FALSE)
save(df_relative_crit, file = paste(base_dir, "/temp/df_relative_crit.RDS",sep=""))