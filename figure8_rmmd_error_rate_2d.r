

# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# Load packages
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(broom)
p_load(gridExtra)
p_load(grid)
p_load(rlang)
p_load(colorspace)
p_load(VGAM)
p_load(boot)
p_load(dplyr)
p_load(cowplot)
p_load(binom)
p_load(VGAM)
p_load(gplots)
p_load(RColorBrewer)
p_load(tidyr)
p_load(docstring)

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
df_results <- quantify_coverage_errors(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma.rds",
                                overwrite=overwrite, mus_a=100, sigmas_a=1)
# Assemble results into square matrix
grid_slopes <- slopes_by_rowcol(df_results$mean_diff_rmmd_rmu, sigmas, mus)
load(file = "temp/debug.RData")




# 2D Error Rate relative x_bar to relative rmu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rxbar_error_rate)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
grid_slopes <- slopes_by_rowcol(df_results$mean_xbar_error_rate, sigmas, mus)
# Plot heatmap
gg <- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
print(gg)
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_1a xbar error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 

# 2D visualization of hypothesized coverage error of rxbar in mu space
#______________________________________________________________________________#
#
# Convert from matrix to data frame
df_err_test <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$pval_rxbar_err_eq_zero > 0.05/(length(sigmas)*length(mus)),
  df_results$pval_rxbar_err_eq_alpha > 0.05/(length(sigmas)*length(mus))))) %>% gather(mu, z, -sigma)
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
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_1b rxbar error test.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)


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
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 0.1), limits=c(0,0.1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_2a rmmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 


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
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_2b mmd error test.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)





# 2D visualization of tested coverage error of rMMD in mu/sigma space
#                                                                              #
#______________________________________________________________________________#
sigmas <- seq(.1, 5, by = .1)
mu_ov_sigmas <- seq (-.5, .5, by=0.01)
n_obs <- 50
set.seed(rand.seed)

# Error rate rMMD
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rmmd_error_rate)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
grid_slopes <- slopes_by_rowcol(df_results$mean_rmmd_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 0.1), limits=c(0,0.1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_2a rmmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 



# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(NULL, sigmas, n_samples, n_obs,
                                "temp/mmd_Error_2D_mu_over_sigma_vs_sigma.rds", mu_ov_sigmas,
                                overwrite = overwrite,mus_a=100, sigmas_a=1)
# load(file = "temp/debug_space.RData")
# Error MMD < mu in mu/sigma space
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$pval_rmmd_err_eq_zero > 0.05/(length(sigmas)*length(mu_ov_sigmas)),
  df_results$pval_rmmd_err_eq_alpha > 0.05/(length(sigmas)*length(mu_ov_sigmas))))) %>% 
  gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu/sigma)) + ylab(expression(sigma)) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_3a mmd error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)

# Error MMD < mu in mu/sigma space
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$pval_rxbar_err_eq_zero > 0.05/(length(sigmas)*length(mu_ov_sigmas)),
  df_results$pval_rxbar_err_eq_alpha > 0.05/(length(sigmas)*length(mu_ov_sigmas))))) %>% 
  gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu/sigma)) + ylab(expression(sigma)) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_3b xbar error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)

# Identify location of coverage error boundaries with mu/sigma space
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
mus <- NULL
sigmas <- seq(1.5, 5, by = .1); 
mu_ov_sigmas <- seq (0.10, 0.40, by=0.001)
source("R/error_2d_utils.R")
df_crit_mu_ov_sigma <- 
  locate_bidir_binary_thresh(ind_var = "rmmd", mus=NULL, sigmas = sigmas, n_obs = n_obs, 
                             temp_name = "mmd_Error_mu_over_sigma_vs_sigma.rds", 
                             mu_ov_sigmas = mu_ov_sigmas, rand.seed = rand.seed,
                             overwrite = overwrite)
df_crit_mu_ov_sigma$merge = paste(df_crit_mu_ov_sigma$er, df_crit_mu_ov_sigma$side)

# Plot of each boundary separate
df_slopes <- df_crit_mu_ov_sigma %>% group_by(er,side) %>% summarize(pearson = cor.test(
  critical_mu_over_sigma, sigma,method = "pearson")$p.value)
df_slopes$adj_pearson <- p.adjust(df_slopes$pearson,"bonferroni")
gg <- ggplot(data=df_crit_mu_ov_sigma,aes(x=sigma, y=critical_mu_over_sigma,
                                          shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  xlab(expression(sigma)) + 
  ylab(expression(mu*phantom(.)*phantom(.)/sigma)) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
# scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  "))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_d mmd boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)


