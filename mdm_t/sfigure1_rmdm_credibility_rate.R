



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

source("R/credibility_rate_toolbox.R")



# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
# Script Parameters
fig_num = "5"
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
xbars_dm <- seq(-2.5, 2.5, by = .05)
sds_dm <- seq(.01, 1, by = .025)
# Assume equal variances
sds_a = sds_dm/sqrt(2/n_obs)
sds_b = sds_a
alpha <- 0.05
n_samples <- 5*250/alpha


# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  process_cred_intervals(xbars_a = 100, sds_a = sds_a, n_a = n_obs, 
                         xbars_b = 100 + xbars_dm, sds_b = sds_b, n_b = n_obs, alpha = alpha,
                         n_samples = n_samples, out_path = paste(fig_path, "/mdm_cred_xbar_vs_s.rds",sep=""),
                         overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                         stat_name = "rmdm", method = "montecarlo")

# 1A: Error rate of MDM < mu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sds_dm, as_tibble(df_results$cred_rate)) %>% gather(mu, z, -sigma)
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
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(1-2*alpha,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.26, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_1a mdm error rate 2D a0p05.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 







# Run simulations calculating error of mdm with mu and sigma swept
alpha=0.025
n_samples <- 5*250/alpha
df_results <- 
  process_cred_intervals(xbars_a = 100, sds_a = sds_a, n_a = n_obs, 
                         xbars_b = 100 + xbars_dm, sds_b = sds_b, n_b = n_obs, alpha = alpha,
                         n_samples = n_samples, out_path = paste(fig_path, "/mdm_cred_xbar_vs_s.rds",sep=""),
                         overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                         stat_name = "rmdm", method = "montecarlo")

# 1A: Error rate of MDM < mu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sds_dm, as_tibble(df_results$cred_rate)) %>% gather(mu, z, -sigma)
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
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(1-2*alpha,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.26, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_1a mdm error rate 2D a0p025.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 






# Run simulations calculating error of mdm with mu and sigma swept
alpha=0.025/2
n_samples <- 5*250/alpha
df_results <- 
  process_cred_intervals(xbars_a = 100, sds_a = sds_a, n_a = n_obs, 
                         xbars_b = 100 + xbars_dm, sds_b = sds_b, n_b = n_obs, alpha = alpha,
                         n_samples = n_samples, out_path = paste(fig_path, "/mdm_cred_xbar_vs_s.rds",sep=""),
                         overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                         stat_name = "rmdm", method = "montecarlo")

# 1A: Error rate of MDM < mu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sds_dm, as_tibble(df_results$cred_rate)) %>% gather(mu, z, -sigma)
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
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(1-2*alpha,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.26, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_1a mdm error rate 2D a0p125.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 





#

n_obs = 6
xbars_dm <- seq(-2.5, 2.5, by = .2)
sds_dm <- seq(.01, 1, by = .05)
# Assume equal variances
sds_a = sds_dm/sqrt(2/n_obs)
sds_b = sds_a

alphas = 0.05/c(1,2.5,5,10,20,50,100, 500)
n_samples = 5*250/alphas

cred_rate = tibble(mean = rep(0, length(alphas)), std = rep(0, length(alphas)), 
                   n = rep(0, length(alphas)), conf.level = 1-alphas)
for (n in seq(1, length(alphas))) {
  
  df_results <- 
    process_cred_intervals(xbars_a = 100, sds_a = sds_a, n_a = n_obs, 
                           xbars_b = 100 + xbars_dm, sds_b = sds_b, n_b = n_obs, alpha = alphas[n],
                           n_samples = n_samples[n], out_path = paste(fig_path, "/mdm_cred_xbar_vs_s.rds",sep=""),
                           overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                           stat_name = "rmdm", method = "montecarlo")
  cred_rate$mean[n] = mean(df_results$cred_rate) 
  cred_rate$std[n] = sd(df_results$cred_rate)
  cred_rate$n[n] = prod(dim(df_results$cred_rate))
}

save(list = ls(all.names = TRUE), file = "temp/sfigrue5_rmdm_credibility.RData",
     envir = environment())
load(file = "temp/sfigrue5_rmdm_credibility.RData")
cred_rate$conf.level = 1-alphas


gg<- ggplot(cred_rate, aes(conf.level, mean)) + 
  geom_point() + geom_line() + xlab("1-alpha") + ylab("Mean Cred. Rate") +
  geom_ribbon(aes(ymin=mean - 1.96*std, ymax=mean - 1.96*std)) +
  coord_cartesian(xlim = c(0.94, 1), ylim = c(0.94, 1)) +
  theme_classic(base_size=8) 
gg
save_plot(paste(fig_path, "/", fig_num, "_1a mdm error rates.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 
ft <- lm(mean~0+conf.level, data = cred_rate)
summary(ft)
confint(ft)




