

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
source("R/equivalence_contests.R")


# Figure parameters
fig_num = "9" 
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/F",fig_num, "/",sep="")
# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_sims = 1e3
n_samples = 1e3
n_obs = 50
rand.seed = 1

parallel_sims = TRUE
include_bf = TRUE
# scale_contest_path = paste("figure/F", fig_num, "/F", fig_num,"_scale_contest_results.csv",sep="")

# Test how each metric responds to changing each characteristic of similarity




# Un-transformed Error
# Contest 1-1) Quantify Error rate with each metric predicting experiment with
# Lower mean difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
set.seed(rand.seed)
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste("F", fig_num, "_1a_esize_contest_mu_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims = 21, rand.seed, 
                                   mus_1a  = rep(1,21), 
                                   sigmas_1a = 0.71, 
                                   mus_1b  = seq(-5,5,0.5), 
                                   sigmas_1b = 0.71,
                                   mus_2a  = rep(1,21), 
                                   sigmas_2a = 0.71, 
                                   mus_2b  = rep(1,21),  
                                   sigmas_2b = 0.71, 
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames)
dfs_1a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path, is_plotted = FALSE)

# Plot stat values over independent variable

df_indvar_sum_1a <- plot_stats_over_var(dfs_1a$df_es, "mus_1b")

df_es <- dfs_1a$df_es
indvar <- "mus_1b"


# var_b = var_b



# Change mu w/ sigma=1, calculate pearson for each stat across mu


# Change sigma at mu=0, calculate pearson for each stat across sigma


# Change sigma at mu=10, calculate pearson for each stat across sigma




# Change rmu w/ sigma=1, calculate pearson for each stat across rmu


# Change rsigma at mu=0, calculate pearson for each stat across rsigma


# Change rsigma at mu=10, calculate pearson for each stat across rsigma


# Summarizing table showing MMD is the only stat with a solid color that is statistically significant from







# CHange sampling with fixed std and mu