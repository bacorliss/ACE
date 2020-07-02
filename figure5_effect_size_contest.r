

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

library(pacman)
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
source("R/effect_size_contests.R")



# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomly chosen
n_sims = 1e3
n_samples = 1e2
n_obs = 50
rand.seed = 0
fig_basename = "F5"
parallel_sims = TRUE
include_bf = TRUE
scale_contest_path = paste("figure/",fig_basename,"_scale_contest_results.csv",sep="")
rscale_contest_path = paste("figure/",fig_basename,"_rscale_contest_results.csv",sep="")


###############################################################################
#
# Un-transformed Error
# Contest 1-1) Quantify Error rate with each metric predicting experiment with
# Lower mean difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste(fig_basename, "_1a_esize_", "contest_mu_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims, 1,5), 
                                   sigmas_1d = rep(15,n_sims),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,1, 5), 
                                   sigmas_2d = rep(15,n_sims),
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)
dfs_1a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                     y_ax_str = "abs(~mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))

# [ Far from zero ]
#
#------------------------------------------------------------------------------
gt_colnames = "is_mud_md2gtmd1"
fig_name = paste(fig_basename, "_1b_esize_contest_mu_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                    sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,1,5), 
                                    sigmas_1d = rep(0.05,n_sims),
                                   mus_2a  = rep(10,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,1,5), 
                                    sigmas_2d = rep(0.05,n_sims),
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)  
dfs_1b <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                 y_ax_str = "abs(~mu[DM]*phantom(.))",
                                 include_bf = include_bf, parallel_sims = parallel_sims,
                                 fig_name = paste(fig_name, ".tiff",sep = ""))


# COntest 2) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means
# [Near from zero]
#
#------------------------------------------------------------------------------
gt_colnames = "is_sigma_md2gtmd1"
fig_name = paste(fig_basename, "_2a_esize_", "contest_sigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(2,n_sims), 
                                   sigmas_1d = runif(n_sims, 5, 30),
                                   
                                   mus_2a  = rep(25,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(2,n_sims), 
                                   sigmas_2d = runif(n_sims, 5, 50),
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames) 
dfs_2a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                     y_ax_str = "sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))


# [Far from zero]
#
#------------------------------------------------------------------------------
gt_colnames = "is_sigma_md2gtmd1"
fig_name = paste(fig_basename, "_2b_esize_contest_sigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(10,n_sims), 
                                   sigmas_1d = runif(n_sims, 5, 20),
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(10,n_sims), 
                                   sigmas_2d = runif(n_sims, 5, 20),
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)
dfs_2b <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                     y_ax_str = "sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))



# Contest 3) Quantify Error rate with each metric predicting experiment with
# Lower STD of difference in means and mean difference in means [Free form]
# [Near from zero]
#
#------------------------------------------------------------------------------
gt_colnames = c("is_mud_md2gtmd1","is_sigma_md2gtmd1")
fig_name = paste(fig_basename, "_3a_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,1,2), 
                                   sigmas_1a = runif(n_sims,2,10), 
                                   mus_1d  = runif(n_sims,0.1,2), 
                                   sigmas_1d = runif(n_sims,2,10),
                                   
                                   mus_2a  = runif(n_sims,50,50),  
                                   sigmas_2a = runif(n_sims,2,10),
                                   mus_2d  = runif(n_sims,0.1,2), 
                                   sigmas_2d = runif(n_sims,2,10),
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames) 
dfs_3a_mu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                     y_ax_str = "abs(~mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_mu.tiff",sep = ""))
dfs_3a_sigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                     y_ax_str = "sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_sigma.tiff",sep = ""))


# [Far from zero]
#
#------------------------------------------------------------------------------
gt_colnames = c("is_mud_md2gtmd1","is_sigma_md2gtmd1")
fig_name = paste(fig_basename, "_3b_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims, 1, 1), 
                                   sigmas_1a = runif(n_sims, 1, 1), 
                                   mus_1d  = runif(n_sims, 5, 10), 
                                   sigmas_1d = runif(n_sims, 5, 15),
                                   
                                   mus_2a  = runif(n_sims, 100, 100),
                                   sigmas_2a = runif(n_sims, 1, 1),
                                   mus_2d  = runif(n_sims, 5, 10), 
                                   sigmas_2d = runif(n_sims, 5, 15),
                                   
                                   switch_group_ab = TRUE,
                                   switch_sign_mean_d = TRUE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames) 
dfs_3b_mu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                     y_ax_str = "abs(~mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_mu.tiff",sep = ""))
dfs_3b_sigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                     y_ax_str = "sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_sigma.tiff",sep = ""))



###############################################################################
#
# Relative Error
#
##############################################################################

# Contest 4) Quantify Error rate with each metric predicting experiment with
# Lower relative mean difference in means
# [Near zero]
#
#------------------------------------------------------------------------------
gt_colnames = "is_rmud_md2gtmd1"
fig_name = paste(fig_basename, "_4a_esize_contest_rmu_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,1,5), 
                                   sigmas_1d = rep(15,n_sims),
                                   
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,5,25), 
                                   sigmas_2d = rep(75,n_sims),
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)  
dfs_4a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                     y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))

# [ Far from zero ]
#
#------------------------------------------------------------------------------
gt_colnames = "is_rmud_md2gtmd1"
fig_name = paste(fig_basename, "_4b_esize_contest_rmu_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = rep(10,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_1d  = runif(n_sims,2,5), 
                                   sigmas_1d = rep(3,n_sims),
                                   
                                   mus_2a  = rep(50,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_2d  = runif(n_sims,10,25), 
                                   sigmas_2d = rep(15,n_sims),
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)  
dfs_4b <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                     y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))


# Contest 5) Quantify Error rate with each metric predicting experiment with
# Lower relative STD of difference in means
# [u_d near zero]
#
#------------------------------------------------------------------------------
gt_colnames = "is_rsigma_md2gtmd1"
fig_name = paste(fig_basename, "_5a_esize_contest_rsigma_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(100,n_sims), 
                                   sigmas_1a = rep(1, n_sims),
                                   mus_1d  = runif(n_sims, 2, 2), 
                                   sigmas_1d = runif(n_sims, 6, 20),
                                   
                                   mus_2a  = rep(500,n_sims), 
                                   sigmas_2a = rep(1, n_sims),
                                   mus_2d  =  runif(n_sims, 10, 10), 
                                   sigmas_2d = runif(n_sims, 30, 100),
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames) 
dfs_5a <- process_esize_simulations(df_init, gt_colname = gt_colnames, 
                                     y_ax_str = "r*sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, ".tiff",sep = ""))
# [ Far from zero ]
#
#------------------------------------------------------------------------------
gt_colnames = "is_rsigma_md2gtmd1"
fig_name = paste(fig_basename, "_5b_esize_contest_rsigma_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(100,n_sims), 
                                    sigmas_1a = rep(1,n_sims),
                                   mus_1d  = rep(20,n_sims), 
                                    sigmas_1d = runif(n_sims, 10, 50),
                                   
                                   mus_2a  = rep(400,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                   mus_2d  = rep(80,n_sims), 
                                    sigmas_2d = runif(n_sims, 40, 200),
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames) 
dfs_5b <- process_esize_simulations(df_init, gt_colname = "is_rsigma_md2gtmd1", 
                      y_ax_str = "r*sigma[DM]",
                      include_bf = include_bf, parallel_sims = parallel_sims,
                      fig_name = paste(fig_name, ".tiff",sep = ""))



# Contest 6) Quantify Error rate with each metric predicting experiment with
# Lower rel. STD, and rel. mean, of difference in means [Free form]
# [Near from zero]
#
#------------------------------------------------------------------------------
gt_colnames = c("is_rmud_md2gtmd1","is_rsigma_md2gtmd1")
fig_name = paste(fig_basename, "_6a_esize_contest_free_near_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,10,10), 
                                   sigmas_1a = runif(n_sims,1,1), 
                                   mus_1d  = runif(n_sims,1,3), 
                                   sigmas_1d = runif(n_sims,8,20),
                                   
                                   mus_2a  = runif(n_sims,40,40),  
                                   sigmas_2a = runif(n_sims,1,1),
                                   mus_2d  = runif(n_sims,4,12),  
                                   sigmas_2d = runif(n_sims,32,80),
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)
dfs_6a_rmu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                     y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rmu.tiff",sep = ""))
dfs_6a_rsigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                     y_ax_str = "r*sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rsigma.tiff",sep = ""))

# [Far from zero]
#
#------------------------------------------------------------------------------
gt_colnames = c("is_rmud_md2gtmd1","is_rsigma_md2gtmd1")
fig_name = paste(fig_basename, "_6b_esize_contest_free_far_zero", sep = "")
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed,
                                   mus_1a  = runif(n_sims,20,20), 
                                   sigmas_1a = runif(n_sims,1,1), 
                                   mus_1d  = runif(n_sims,4,8), 
                                   sigmas_1d = runif(n_sims,4,12),
                                   
                                   mus_2a  = runif(n_sims,80,80),  
                                   sigmas_2a = runif(n_sims,1,1),
                                   mus_2d  = runif(n_sims,16,32),  
                                   sigmas_2d = runif(n_sims,16,48),
                                   
                                   switch_sign_mean_d = TRUE,
                                   switch_sign_mean_ab = FALSE,
                                   switch_group_ab = FALSE,
                                   fig_name = paste(fig_name, ".tiff",sep = ""),
                                   gt_colnames=gt_colnames)
dfs_6b_rmu <- process_esize_simulations(df_init, gt_colname = gt_colnames[1], 
                                     y_ax_str = "abs(~r*mu[DM]*phantom(.))",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rmu.tiff",sep = ""))
dfs_6b_rsigma <- process_esize_simulations(df_init, gt_colname = gt_colnames[2], 
                                     y_ax_str = "r*sigma[DM]",
                                     include_bf = include_bf, parallel_sims = parallel_sims,
                                     fig_name = paste(fig_name, "_rsigma.tiff",sep = ""))




# Save and load effect size contest data
save.image(file=paste("figure/", fig_basename, "_effect_size_contest_results.png"))
load(paste("figure/", fig_basename, "_effect_size_contest_results.png"))





# OVerall contest heatmaps
#------------------------------------------------------------------------

# Re-make heatmap with rectangles based on the selection
my_palette <- colorRampPalette(c(rgb(255, 0, 0,maxColorValue = 255),"white",
                                 rgb(47, 117, 181,maxColorValue = 255)))(n = 299)
col_breaks = c(seq(-1, -.1, length=100), seq(-.09, 0.09, length=100), 
               seq(0.1, 1.0,length=100))
# Function for making selection rectangles around selection cells
makeRects <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[2])[cells,]
  xl=coords[,2]-0.49; yb=coords[,1]-0.49; xr=coords[,2]+0.49; yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}
make2RectGroups <- function(cells1,lwd1, cells2, lwd2){
  makeRects(cells1,lwd1)
  makeRects(cells2,lwd2)
}




# Export summary stats for scale data
#-------------------------------------------------------------------------------
dfs_scale <- c(ls(pat="dfs_[1-3]a"),ls(pat="dfs_[1-3]b"))
scale_norm_ind = rep(c(1,3,1,3),2)

# Extract means for each group and subtract from 0.5 (random)
scale_means_from_0.5 <- 0.5 - sapply(dfs_scale, function(x) get(x)$df_plotted$mean)
scale_scores <- t(t(scale_means_from_0.5)/
                    scale_means_from_0.5[cbind(scale_norm_ind,seq_along(scale_norm_ind))])
# Export csv
rownames(scale_scores) <- effect_size_dict[[4]]
# Get statistical significance
scale_scores_sig <- !sapply(dfs_scale, function(x) get(x)$df_plotted$is_mean_0.5) 
scale_score_norm <- sapply(scale_norm_ind, function(ind,len) 
  ifelse(1:len == ind, TRUE,FALSE), length(effect_size_dict[[4]]))

png(paste("figure/", fig_basename, "es_contest scale.png"),    
    width = 5*300, height = 2.5*300, res = 300, pointsize = 8)  
heatmap.2(scale_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {make2RectGroups(scale_scores_sig,1,scale_score_norm,3)}, 
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(effect_size_dict[[4]], function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(scale_scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(scale_scores)[1]),
          notecol="black",notecex=1, lwid=c(0.01,5),lhei=c(0.01,5),margins =c(.1,.1))
dev.off()

# Export summary stats for relative scale data
#-------------------------------------------------------------------------------
# Export summary stats for relative scale data
dfs_rscale <- c(ls(pat="dfs_[4-6]a"),ls(pat="dfs_[4-6]b"))
rscale_norm_ind = rep(c(2,4,2,4),2)
# Extract means for each group and subtract from 0.5 (random)
rscale_means_from_0.5 <- 0.5 - sapply(dfs_rscale, function(x) get(x)$df_plotted$mean)
rscale_scores <- t(t(rscale_means_from_0.5)/
                     rscale_means_from_0.5[cbind(rscale_norm_ind,seq_along(rscale_norm_ind))])
# Export csv
rownames(rscale_scores) <- effect_size_dict[[2]]
# Get statistical significance
rscale_scores_sig <- !sapply(dfs_rscale, function(x) get(x)$df_plotted$is_mean_0.5) 
rscale_score_norm <- sapply(rscale_norm_ind, function(ind,len) 
  ifelse(1:len == ind, TRUE,FALSE), length(effect_size_dict[[4]]))

png(paste("figure/", fig_basename, "es_contest relative scale.png"),    
    width = 5*300, height = 2.5*300, res = 300, pointsize = 8)  
heatmap.2(rscale_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {make2RectGroups(rscale_scores_sig,1,rscale_score_norm,3)}, 
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(effect_size_dict[[4]], function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(rscale_scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(rscale_scores)[1]),
          notecol="black",notecex=1, lwid=c(0.01,5),lhei=c(0.01,5),margins =c(.1,.1))
dev.off()



# Contest 7) which metric is best at discerning experiments below a threshold 
# difference in means
#
#------------------------------------------------------------------------------




# Contest 8) which metric is best at discerning experiments below a threshold 
# relative difference in means
#
#------------------------------------------------------------------------------




