

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
fig_num = "1" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")


# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen

n_samples = 1e3
n_obs = 50
rand.seed = 1
gt_colnames = "is_mudm_2gt1"
parallel_sims = TRUE
include_bf = TRUE
# scale_contest_path = paste("figure/F", fig_num, "/F", fig_num,"_scale_contest_results.csv",sep="")

# Test how each metric responds to changing each characteristic of similarity





# Un-scaled Space
#
#-------------------------------------------------------------------------------

# Mu
# Pearson rho of abs(mu_d1) versus abs(mean of each stat)
# Sweep mu_1d, sigma_d = 1
#------------------------------------------------------------------------------
# Fixed mu_d, but as it increases, rmu_d decreases
set.seed(rand.seed)
mus_d_vect = seq(4.9,0.1,-0.2)
mus_a_vect = mus_d_vect*2
mus_b_vect = mus_d_vect + mus_a_vect; n_sims = length(mus_b_vect)
sigmas_ab_vect=1e-2

gt_colnames = "is_mudm_2gt1"
fig_name = paste("F", fig_num, "_1a_stat_correlation_mu_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples = n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_a_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_b_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(0,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)

df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = parallel_sims,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mu_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "mu_1dm", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_agreement = 1,
                           stats_basenames = attr(df_esize$df_es,"varnames"),
                           stats_labels = attr(df_esize$df_es,"varnames_pretty"))






# Sigma
# Pearson rho of sigma versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigmas_ab_vect = seq(10,1,-0.25); n_sims = length(sigmas_ab_vect)
mus_a_vect = sigmas_ab_vect
mus_b_vect = mus_a_vect;
  

gt_colnames = "is_mudm_2gt1"
fig_name = paste("F", fig_num, "_1b_stat_correlation_sigma_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_a_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_b_vect, 
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "sigma[pool]",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_sigma_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "sigma_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_agreement = 1,
                           stats_basenames = attr(df_esize$df_es,"varnames"),
                           stats_labels = attr(df_esize$df_es,"varnames_pretty"))


# Sample Size:   Pearson rho of sigma versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
# TODO not working properly, I want to make rsigmaDM a constant, can't make it work
#------------------------------------------------------------------------------
set.seed(rand.seed)
n_1ab_vect = seq(5, 50,  2); n_sims = length(n_1ab_vect)
mus_ab_vect = 1
sigmas_ab_vect = 1

gt_colnames = "is_mudm_2gt1"
fig_name = paste("F", fig_num, "_1c_stat_correlation_df_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_ab_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_ab_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(0,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_1ab_vect, n_1b = n_1ab_vect, n_2a = 30, n_2b = 30,
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "df[pool]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_df_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "df_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, dir_to_agreement = -1,
                           stats_basenames = attr(df_esize$df_es,"varnames"),
                           stats_labels = attr(df_esize$df_es,"varnames_pretty"))







# Relative Scale

#------------------------------------------------------------------------------

# Relative Mean
#------------------------------------------------------------------------------
# Decreasing Rmu is more similar

# Increase mu a
# Same d
# Increase std for a and b to match mu a
set.seed(rand.seed)
mus_a_vect =  seq(10,20,0.5); n_sims = length(mus_a_vect) 
mus_b_vect =  mus_a_vect+1
sigmas_ab_vect = mus_a_vect + 0.5 * (mus_b_vect - mus_a_vect)



gt_colnames = "is_mudm_2gt1"
fig_name = paste("F", fig_num, "_1d_stat_correlation_rmu_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed,
                                      mus_1a  = mus_a_vect,
                                      sigmas_1a = sigmas_ab_vect,
                                      mus_1b  = mus_b_vect,
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(0,n_sims),
                                      sigmas_2a = 1,
                                      mus_2b  = rep(0,n_sims),
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)

df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames,
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_rmu_pearson <-
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "rmu_1dm",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_agreement = 1,
                           stats_basenames = attr(df_esize$df_es,"varnames"),
                           stats_labels = attr(df_esize$df_es,"varnames_pretty"))



# Relative sigma
#------------------------------------------------------------------------------
set.seed(rand.seed)
mus_b_vect =  seq(10,20,0.5); n_sims = length(mus_b_vect) 
mus_a_vect = mus_b_vect
sigmas_ab_vect = 1

gt_colnames = "is_mudm_2gt1"
fig_name = paste("F", fig_num, "_1e_stat_correlation_rsigma_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_a_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_b_vect, 
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(0,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(0,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      switch_group_ab = FALSE,
                                      switch_sign_mean_d = FALSE,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_rsigma_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "rsigma_1dm", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, dir_to_agreement = 1,
                           stats_basenames = attr(df_esize$df_es,"varnames"),
                           stats_labels = attr(df_esize$df_es,"varnames_pretty"))


# Sample Size:   Pearson rho of sigma versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
n_1ab_vect = seq(16, 100, 2); n_sims = length(n_1ab_vect)
mus_a_vect = 1
mus_b_vect = 1
sigmas_ab_vect = 1

fig_name = paste("F", fig_num, "_1f_stat_correlation_rdf_sweep", sep = "")
df_init <- generateExperiment_Data(n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_a_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_b_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(0,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(0,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_1ab_vect, n_1b = n_1ab_vect, n_2a = 30, n_2b = 30,
                                   switch_group_ab = FALSE,
                                   switch_sign_mean_d = FALSE,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "sigma[pool]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_rdf_pearson <- 
  lineplot_indvar_vs_stats(df = df_esize$df_es, indvar = "df_1d", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_agreement = -1,
                           stats_basenames = attr(df_esize$df_es,"varnames"),
                           stats_labels = attr(df_esize$df_es,"varnames_pretty"))










# Summary Heatmaps
#------------------------------------------------------------------------

# Re-make heatmap with rectangles based on the selection
my_palette <- colorRampPalette(c(rgb(255, 0, 0,maxColorValue = 255),"white",
                                 rgb(47, 117, 181,maxColorValue = 255)))(n = 299)
col_breaks = c(seq(-1, -.1, length=100), seq(-.09, 0.09, length=100), 
               seq(0.1, 1.0,length=100))
# Function for making selection rectangles around selection cells
makeRects <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[1])[cells,]
  xl=coords[,2]-0.49; yb=coords[,1]-0.49; xr=coords[,2]+0.49; yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}
make2RectGroups <- function(cells1,lwd1){
  makeRects(cells1,lwd1)
}



# Unscaled Heatmap Summary
#-------------------------------------------------------------------------------
scores = cbind(df_mu_pearson$pearson_rho, df_sigma_pearson$pearson_rho,
               df_df_pearson$pearson_rho)
scores_sig = cbind(df_mu_pearson$is_pearson_rho_sig, df_sigma_pearson$is_pearson_rho_sig,
                   df_df_pearson$is_pearson_rho_sig)
png(paste(base_dir, "/figure/F", fig_num, "/F", fig_num, "_pearson_unscaled_units.png",sep=""),    
    width = 1.75*300, height = 2*300, res = 300, pointsize = 8)  
heatmap.2(scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {make2RectGroups(scores_sig,2)},
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(200,200),sepcolor="white",
          labRow =  sapply(effect_size_dict[[4]], function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(scores)[1]),
          notecol="black",notecex=1, lwid=c(0.001,5),lhei=c(0.001,5),margins =c(0,0))
dev.off()






# Relative Heatmap Summary
#-------------------------------------------------------------------------------
scores = cbind(df_rmu_pearson$pearson_rho, df_rsigma_pearson$pearson_rho,
               df_rdf_pearson$pearson_rho)
scores_sig = cbind(df_rmu_pearson$is_pearson_rho_sig, df_rsigma_pearson$is_pearson_rho_sig,
                   df_rdf_pearson$is_pearson_rho_sig)
png(paste(base_dir, "/figure/F", fig_num, "/F", fig_num, "_pearson_relative_scale_units.tif",sep=""),    
    width = 1.75*300, height = 2*300, res = 300, pointsize = 8)  
heatmap.2(scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {make2RectGroups(scores_sig,2)},
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(effect_size_dict[[4]], function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(scores)[1]),
          notecol="black",notecex=1, lwid=c(0.001,5),lhei=c(0.001,5),margins =c(0,0))
dev.off()
