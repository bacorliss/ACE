

# p.adjust(p, method = p.adjust.methods, n = length(p))


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
# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
# library(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
# library(TOSTER)
library(docstring)

# Script Parameters
fig_num = "2"
dir.create(file.path(getwd(), paste("figure/SF",fig_num,sep="")), showWarnings = FALSE)
n_samples <- 1e3
rand.seed <- 0
overwrite = FALSE

# Helper Functions
source("R/mmd.R")
source("R/error_2d_utils.R")






# 2D visualization of mmd difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
n_obs <- 50
# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma.rds",
                                overwrite = overwrite) 

# Difference MMD from Mu
#------------------------------------------------------------------------------
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_diff_mmd_mu)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
# Assemble results into square matrix
grid_slopes <- rowcol_slopes(df_results$mean_diff_mmd_mu, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  # annotate("text", x=max(mus), y=sigmas+.02, label = df_slopes$row_sig_labels, size=2,vjust=1) +
  # annotate("text", x=mus, y=max(sigmas)+.02, label = df_slopes$col_sig_labels, size=2,vjust=1) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(0,2)) +
  # coord_cartesian(xlim = c(min(mus)+0.1, max(mus)+0.1), # This focuses the x-axis on the range of interest
  #                 clip = 'off') +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a1 mmd diff.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 
# Plot slope of heatmap by column
gg <- ggplot(data = grid_slopes$df_col, aes(x = col_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_col$slope_95L, ymax = grid_slopes$df_col$slope_95U), alpha=0.2) +
  xlab(expression(mu)) + ylab(expression(paste("Slope By Column  (Diff ~ ",sigma,")"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a1 mmd diff_vs_mus.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 
# Plot slope of heatmap by row
gg <- ggplot(data = grid_slopes$df_row, aes(x = row_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_row$slope_95L, ymax = grid_slopes$df_row$slope_95U), alpha=0.2) +
  xlab(expression(sigma)) + ylab(expression(paste("Slope By Row  (Diff ~ ",abs(~mu*phantom(.))*phantom(.)," )"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a1 mmd diff_vs_sigmas.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 





# Relative Difference MMD normalized by sigma
#------------------------------------------------------------------------------
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rdiff_mmd_mu_ov_sigma)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
grid_slopes <- rowcol_slopes(df_results$mean_rdiff_mmd_mu_ov_sigma, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  # annotate("text", x=max(mus), y=sigmas+.02, label = df_slopes$row_sig_labels, size=2,vjust=1) +
  # annotate("text", x=mus, y=max(sigmas)+.02, label = df_slopes$col_sig_labels, size=2,vjust=1) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a2 mmd rdiff_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 
# Plot slope of heatmap by column
gg <- ggplot(data = grid_slopes$df_col, aes(x = col_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_col$slope_95L, ymax = grid_slopes$df_col$slope_95U), alpha=0.2) +
  xlab(expression(mu)) + ylab(expression(paste("Slope By Column  (RDiff ~ ",sigma,")"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a2 mmd rdiff_over_sigma_col.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 
# Plot slope of heatmap by row
gg <- ggplot(data = grid_slopes$df_row, aes(x = row_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_row$slope_95L, ymax = grid_slopes$df_row$slope_95U), alpha=0.2) +
  xlab(expression(sigma)) + ylab(expression(paste("Slope By Row  (RDiff ~ ",abs(~mu*phantom(.))*phantom(.)," )"))) +
  theme_classic(base_size=8)
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a2 mmd rdiff_over_sigma_row.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 






  
# Relative Difference MMD normalized by sigma
#------------------------------------------------------------------------------
# mus <- seq(-2.5, 2.5, by = .1); mus<- mus[mus!=0]
# sigmas <- seq(.1, 5, by = .1)
# n_obs <- 50
# # Run simulations calculating error of mmd with mu and sigma swept
# df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma_no_zero.rds",
#                                 overwrite = TRUE) 
source("R/error_2d_utils.R")
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rdiff_mmd_mu_ov_mu)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
grid_slopes <- rowcol_slopes(df_results$mean_rdiff_mmd_mu_ov_mu, sigmas, mus)
# Plot heatmap
gg <- ggplot(df, aes(mu, sigma, fill = z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  # annotate("text", x=max(mus), y=sigmas+.02, label = df_slopes$row_sig_labels, size=2,vjust=1) +
  # annotate("text", x=mus, y=max(sigmas)+.02, label = df_slopes$col_sig_labels, size=2,vjust=1) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a2 mmd rdiff_over_mu.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 
# Plot slope of heatmap by column
gg <- ggplot(data = grid_slopes$df_col, aes(x = col_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_col$slope_95L, ymax = grid_slopes$df_col$slope_95U), alpha=0.2) +
  xlab(expression(mu)) + ylab(expression(paste("Slope By Column  (RDiff ~ ",sigma,")"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a2 mmd rdiff_over_mu_col.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 
# Plot slope of heatmap by row
gg <- ggplot(data = grid_slopes$df_row, aes(x = row_vals, y = slope)) +
  geom_line(size=0.5) + geom_point(size=0.4) + 
  geom_ribbon(aes(ymin = grid_slopes$df_row$slope_95L, ymax = grid_slopes$df_row$slope_95U), alpha=0.2) +
  xlab(expression(sigma)) + ylab(expression(paste("Slope By Row  (RDiff ~ ",abs(~mu*phantom(.))*phantom(.)," )"))) +
  theme_classic(base_size=8)
gg
save_plot(paste("figure/SF", fig_num, "/SF", fig_num, "_a2 mmd rdiff_over_mu_row.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 

