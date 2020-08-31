



# Load packages
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)
library(rlang)
library(colorspace)
library(VGAM)
library(boot)
library(dplyr)
library(cowplot)
library(binom)
library(VGAM)
library(gplots)
library(RColorBrewer)
library(tidyr)
library(docstring)

# Script Parameters
fig_num = "8"
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
rand.seed <- 0
overwrite <- TRUE

# Helper Functions
source("R/mmd.R")
source("R/error_2d_utils.R")


# 2D visualization of mmd difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
n_obs <- 50
n_samples <- 1e3

# debug(stats_param_sweep)
# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma.rds",
                                overwrite=overwrite, mus_a=100, sigmas_a=1)
# Assemble results into square matrix
grid_slopes <- slopes_by_rowcol(df_results$mean_diff_rmmd_rmu, sigmas, mus)


# Error rate of MMD < Mu
#------------------------------------------------------------------------------
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rmmd_error_rate)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
grid_slopes <- slopes_by_rowcol(df_results$mean_rmmd_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a1 rmmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 
# Plot slope of heatmap by column
gg <- ggplot(data = grid_slopes$df_col, aes(x = col_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_col$slope_95L, ymax = grid_slopes$df_col$slope_95U), alpha=0.2) +
  xlab(expression(mu)) + ylab(expression(paste("Slope By Row  (ER ~ ",sigma,")"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a2 rmmd error rate 2D_vs_mus.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 
# Plot slope of heatmap by row
gg <- ggplot(data = grid_slopes$df_row, aes(x = row_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_row$slope_95L, ymax = grid_slopes$df_row$slope_95U), alpha=0.2) +
  xlab(expression(sigma)) + ylab(expression(paste("Slope By Row  (ER ~ ",abs(~mu*phantom(.))," )"))) + 
  theme_classic(base_size=8)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a3 rmmd error rate 2D_vs_sigmas.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 



# 2D visualization of hypothesized coverage error of MMD in mu space
#______________________________________________________________________________#
#
# Convert from matrix to data frame
df_err_test <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$pval_rmmd_err_eq_zero > 0.05/(length(sigmas)*length(mus)),
  df_results$pval_rmmd_err_eq_alpha > 0.05/(length(sigmas)*length(mus))))) %>% gather(mu, z, -sigma)
df_err_test$mu <- as.numeric(df_err_test$mu)
df_err_test$sigma <- as.numeric(df_err_test$sigma)
df_err_test$z <- factor(df_err_test$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df_err_test, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_b4 mmd error test.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)



# 2D Error Rate x_bar to mu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rxbar_error_rate)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
grid_slopes <- slopes_by_rowcol(df_results$mean_xbar_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a1 rmmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 
# Plot slope of heatmap by column
gg <- ggplot(data = grid_slopes$df_col, aes(x = col_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_col$slope_95L, ymax = grid_slopes$df_col$slope_95U), alpha=0.2) +
  xlab(expression(mu)) + ylab(expression(paste("Slope By Row  (ER ~ ",sigma,")"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a2 rmmd error rate 2D_vs_mus.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 
# Plot slope of heatmap by row
gg <- ggplot(data = grid_slopes$df_row, aes(x = row_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_row$slope_95L, ymax = grid_slopes$df_row$slope_95U), alpha=0.2) +
  xlab(expression(sigma)) + ylab(expression(paste("Slope By Row  (ER ~ ",abs(~mu*phantom(.))," )"))) + 
  theme_classic(base_size=8)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a3 rmmd error rate 2D_vs_sigmas.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 




