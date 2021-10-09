



#' QUantifies coverage error of the mdm (defined as how often abs(mu)>mdm_95 
#' from repeated samples)
#' Results computed in a grid with mu and sigma swept on each exis respectively

# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# Load packages
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(dplyr)
p_load(cowplot)

source("R/credbility_testing.R")



# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
# Script Parameters
fig_num = "4"
fig_path = file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE,recursive = TRUE)
n_samples <- 1e3
rand.seed <- 0
overwrite <- TRUE
is_parallel_proc <- TRUE




# Row 1: 2D visualization of mdm difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#

n_obs = 6
xbars_dm <- seq(-2.5, 2.5, by = .1)
sds_dm <- seq(.01, 1, by = .05)

# Spread sigma_dm across sigma_a and sigma_b equally
sds_a = sds_dm/sqrt(2/n_obs)
sds_b = sds_a
n_samples <- 1e3


# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  process_cred_intervals(xbars_a = 100, sds_a = sds_a, n_a = n_obs, 
                        xbars_b = 100 + xbars_dm, sds_b = sds_b, n_b = n_obs, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/mdm_cred_xbar_vs_s.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                        included_stats = c("mdm"))




# 1A: Error rate of MDM < mu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sds_dm, as_tibble(df_results$credibility_rate_mu_dm_lt_mdm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
# grid_slopes <- slopes_by_rowcol(df_results$mean_mdm_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c(rep("blue",17),"white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(0,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_1a mdm error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 
