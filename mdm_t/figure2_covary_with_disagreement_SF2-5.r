


#' Sample a series of population parameter configurations of a control and experiment group
#' where one parameter of agreement is swept towards increasing disagreement. 
#' Candidate statistics means are calculated based on repeated samples.
#' Correlation between mean value of each statistic and the agreement parameters 
#' are calculated and visualized in a heat map table. 

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
# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
fig_num = "1" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
# Simulation parameters
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_samples = 1e3
n_obs = 6
rand.seed = 1
gt_colnames = "is_mudm_1ldt2"
parallel_sims = TRUE
include_bf = TRUE



# Raw Space
#
#-------------------------------------------------------------------------------

# Unscaled Mu: spearman rho of abs(mu_d1) versus abs(mean of each stat)
#------------------------------------------------------------------------------
# Fixed mu_d, but as it increases, rmu_d decreases
set.seed(rand.seed)
mus_d_vect = seq(5.00, .1,-0.25)
mus_a_vect = mus_d_vect
mus_b_vect = mus_d_vect + mus_a_vect; n_sims = length(mus_b_vect)
sigmas_ab_vect = 1e-4

gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_raw_mu", sep = "")
df_init <- generate_population_configs(n_samples = n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_a_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_b_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(10,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(10,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                    y_ax_str = "abs(~mu[DM]*phantom(.))",
                                    include_bf = include_bf, parallel_sims = TRUE,
                                    fig_name = paste(fig_name, ".tiff",sep = ""),
                                    fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_mu_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "mu_1dm",  indvar_pretty = "mu[DM]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_better = 1)


# Unscaled Sigma: spearman rho of sigma versus abs(mean of each stat)
#------------------------------------------------------------------------------
set.seed(rand.seed)
sigmas_ab_vect = seq(10,1,-0.25); n_sims = length(sigmas_ab_vect)
mus_a_vect = sigmas_ab_vect*10
mus_b_vect = mus_a_vect;
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_raw_sigma", sep = "")
df_init <- generate_population_configs(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_a_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_b_vect, 
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(10,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(10,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "sigma[D]",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_sigma_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "sigma_1d", indvar_pretty = "sigma[D]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_better = 1)


# Unscaled Sample Size:   spearman rho of sigma versus abs(mean of each stat)
# Sweep sigma_1d, mu_1d = 10
#------------------------------------------------------------------------------
set.seed(rand.seed)
n_1ab_vect = seq(5, 50,  2); n_sims = length(n_1ab_vect)
mus_ab_vect = 10
sigmas_ab_vect = .5
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_raw_df", sep = "")
df_init <- generate_population_configs(n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_ab_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_ab_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(10,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(10,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_1ab_vect, n_1b = n_1ab_vect, n_2a = 30, n_2b = 30,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "df[D]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_df_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "df_1d",indvar_pretty = "df[D]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, dir_to_better = -1)

# Unscaled Alpha:   increasing alpha increases agreement
#------------------------------------------------------------------------------
set.seed(rand.seed)
alpha_1 = 0.05/seq(1, 20,1)
alpha_2 = 0.05/seq(1, 20,1)
n_sims = length(alpha_1)
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_raw_alpha", sep = "")
df_init <- generate_population_configs(n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = 20, 
                                   sigmas_1a = 1, 
                                   mus_1b  = 20, 
                                   sigmas_1b = 1,
                                   mus_2a  = 10, 
                                   sigmas_2a = 1, 
                                   mus_2b  = 10,  
                                   sigmas_2b = 1,
                                   n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                   alpha_1 = alpha_1,
                                   alpha_2 = alpha_2,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "Alpha[DM]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_alpha_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "alpha_1",  indvar_pretty = "alpha[DM]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_better = -1)






# Relative Scale

#------------------------------------------------------------------------------

# Relative Mean:  decreasing rmu has higher agreement
#------------------------------------------------------------------------------
set.seed(rand.seed)
mus_a_vect =  seq(10,20*10,0.2*10); n_sims = length(mus_a_vect) 
mus_b_vect =  mus_a_vect+10
sigmas_ab_vect = 0.1*mus_a_vect 
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_rel_mu", sep = "")
df_init <- generate_population_configs(n_samples, n_sims = n_sims, rand.seed,
                                      mus_1a  = mus_a_vect,
                                      sigmas_1a = sigmas_ab_vect,
                                      mus_1b  = mus_b_vect,
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(10,n_sims),
                                      sigmas_2a = 1,
                                      mus_2b  = rep(10,n_sims),
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)

df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames,
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_rmu_spearman <-
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "rmu_1dm", indvar_pretty = "r*mu[DM]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_better = 1)


# Relative sigma: decreasing rsigma has higher agreement
#------------------------------------------------------------------------------
set.seed(rand.seed)
mus_b_vect =  seq(10,20,0.5); n_sims = length(mus_b_vect) 
mus_a_vect = mus_b_vect
sigmas_ab_vect = 1
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_rel_rsigma", sep = "")
df_init <- generate_population_configs(n_samples, n_sims = n_sims, rand.seed, 
                                      mus_1a  = mus_a_vect, 
                                      sigmas_1a = sigmas_ab_vect, 
                                      mus_1b  = mus_b_vect, 
                                      sigmas_1b = sigmas_ab_vect,
                                      mus_2a  = rep(10,n_sims), 
                                      sigmas_2a = 1, 
                                      mus_2b  = rep(10,n_sims),  
                                      sigmas_2b = 1,
                                      n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                      fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                      gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                         y_ax_str = "abs(~mu[DM]*phantom(.))",
                                         include_bf = include_bf, parallel_sims = parallel_sims,
                                         fig_name = paste(fig_name, ".tiff",sep = ""),
                                         fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_rsigma_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "rsigma_1d", indvar_pretty = "r*sigma[DM]", 
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path, dir_to_better = 1)


# Relative Sample Size:   increasing df has higher agreement
#------------------------------------------------------------------------------
set.seed(rand.seed)
# n_1ab_vect = seq(16, 100, 2); n_sims = length(n_1ab_vect)
n_1ab_vect = seq(44, 100, 2); n_sims = length(n_1ab_vect)
mus_a_vect = 10
mus_b_vect = 10
sigmas_ab_vect = 1
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_rel_df", sep = "")
df_init <- generate_population_configs(n_samples, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = mus_a_vect, 
                                   sigmas_1a = sigmas_ab_vect, 
                                   mus_1b  = mus_b_vect, 
                                   sigmas_1b = sigmas_ab_vect,
                                   mus_2a  = rep(10,n_sims), 
                                   sigmas_2a = 1, 
                                   mus_2b  = rep(10,n_sims),  
                                   sigmas_2b = 1,
                                   n_1a = n_1ab_vect, n_1b = n_1ab_vect, n_2a = 30, n_2b = 30,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "sigma[D]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_rdf_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "df_1d", indvar_pretty = "df[D]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_better = -1)
# 

# Relative Alpha:   increasing alpha increases agreement
#------------------------------------------------------------------------------
set.seed(rand.seed)
alpha_1 = 0.05/seq(1, 10,0.5)
n_sims = length(alpha_1)
gt_colnames = "is_mudm_1ldt2"
fig_name = paste("F", fig_num, "_stat_correlation_rel_alpha", sep = "")
df_init <- generate_population_configs(n_samples = 1e2, n_sims = n_sims, rand.seed = rand.seed, 
                                   mus_1a  = 50, 
                                   sigmas_1a = 1, 
                                   mus_1b  = 10, 
                                   sigmas_1b = 1,
                                   mus_2a  = 50, 
                                   sigmas_2a = 1, 
                                   mus_2b  = 10,  
                                   sigmas_2b = 1,
                                   n_1a = n_obs, n_1b = n_obs, n_2a = n_obs, n_2b = n_obs,
                                   alpha_1 = alpha_1,
                                   alpha_2 = alpha_1,
                                   fig_name = paste(fig_name, ".tiff", sep = ""), fig_path = fig_path,
                                   gt_colnames = gt_colnames, is_plotted = FALSE)
df_esize <- process_agreement_contest(df_init, gt_colname = gt_colnames, 
                                      y_ax_str = "Alpha[DM]",
                                      include_bf = include_bf, parallel_sims = parallel_sims,
                                      fig_name = paste(fig_name, ".tiff",sep = ""),
                                      fig_path = fig_path, is_plotted = FALSE)
# Plot stat values over independent variable
df_ralpha_spearman <- 
  plot_stats_covary_indvar(df = df_esize$df_es, indvar = "alpha_1", indvar_pretty = "alpha[DM]",
                           fig_name = paste(fig_name, ".tiff",sep = ""),
                           fig_path = fig_path,  dir_to_better = -1)






# Summary Heatmaps
#------------------------------------------------------------------------

# Re-make heatmap with rectangles based on the selection
my_palette <- colorRampPalette(c(rgb(47, 117, 181,maxColorValue = 255),"white",
                               rgb(255, 0, 0,maxColorValue = 255)))(n = 299)
col_breaks = c(seq(-1, -.1, length=100), seq(-.09, 0.09, length=100), 
               seq(0.1, 1.0,length=100))
# Function for making selection rectangles around selection cells
makeRects <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[1])[cells,]
  xl=coords[,2]-0.45; yb=coords[,1]-0.45; xr=coords[,2]+0.45; yt=coords[,1]+0.45
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}
make2RectGroups <- function(cells1,lwd1){
  makeRects(cells1,lwd1)
}

add_underline <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[2])[cells,]
  xl=coords[,2]-0.49; yb=coords[,1]-0.49; xr=coords[,2]+0.49; yt=coords[,1]+0.49
  segments(xl+.05, yb+.16, xr-.05, yb+.16, col = "black", lty = "solid", lwd = lwd)
}


# Unscaled Heatmap Summary
#-------------------------------------------------------------------------------
scores = cbind(df_mu_spearman$spearman_rho, df_sigma_spearman$spearman_rho,
               df_df_spearman$spearman_rho, df_alpha_spearman$spearman_rho)
scores_sig = cbind(df_mu_spearman$is_spearman_rho_sig, df_sigma_spearman$is_spearman_rho_sig,
                   df_df_spearman$is_spearman_rho_sig, df_alpha_spearman$is_spearman_rho_sig)
# Zero color to white for fields that are not statistically significant
zeroed_scores = scores
zeroed_scores[!scores_sig] <- 0

png(paste(base_dir, "/figure/F", fig_num, "/F", fig_num, "_spearman_unscaled_units.png",sep=""),    
    width = 1.75*300, height = 2*300, res = 300, pointsize = 8)  
heatmap.2(zeroed_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {add_underline(scores_sig,1.5)},
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(200,200),sepcolor="white",
          labRow =  sapply(attr(df_esize,"varnames_pretty"), function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(scores)[1]),
          notecol="black",notecex=1, lwid=c(0.001,5),lhei=c(0.001,5),margins =c(0,0))
dev.off()




# Relative Heatmap Summary
#-------------------------------------------------------------------------------
scores = cbind(df_rmu_spearman$spearman_rho, df_rsigma_spearman$spearman_rho,
               df_rdf_spearman$spearman_rho, df_ralpha_spearman$spearman_rho)
scores_sig = cbind(df_rmu_spearman$is_spearman_rho_sig, df_rsigma_spearman$is_spearman_rho_sig,
                   df_rdf_spearman$is_spearman_rho_sig, df_ralpha_spearman$is_spearman_rho_sig)
# Zero color to white for fields that are not statistically significant
zeroed_scores = scores
zeroed_scores[!scores_sig] <- 0

png(paste(base_dir, "/figure/F", fig_num, "/F", fig_num, "_spearman_relative_scale_units.tif",sep=""),    
    width = 1.75*300, height = 2*300, res = 300, pointsize = 8)  
heatmap.2(zeroed_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {add_underline(scores_sig,1.5);},
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(attr(df_esize,"varnames_pretty"), function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(scores)[1]),
          notecol="black",notecex=1, lwid=c(0.001,5),lhei=c(0.001,5),margins =c(0,0))
dev.off()
