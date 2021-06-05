

#' Quantifies coverage error of the relative mdm (defined as how often abs(rmu)>rmdm_95 
#' from repeated samples)
#' Results computed in a grid with mu and sigma swept on each exis respectively

# Load required packages
#-------------------------------------------------------------------------------
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
p_load(stringr)
# User toolboxes
source("R/coverage_error_toolbox.R")

# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_z"
# Script Parameters
fig_num = "5"
fig_path = file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE, recursive=TRUE)
rand.seed <- 0
overwrite <- TRUE




source("R/coverage_error_toolbox.R")
# Investigation 1: Initial coverage error of rmdm in mu space
#
#_______________________________________________________________________________
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
mus_dm <- seq(-2, 2, by = .1)
sigmas_dm <- seq(.01, .5, by = .025)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_ab = sigmas_dm/sqrt(2/n_obs)
n_samples <- 1e3
mu_ov_sigmas = NULL
# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  quant_coverage_errors(mus_a = mus_a, sigmas_a = sigmas_ab, n_a = n_obs, 
                        mus_b = mus_a + mus_dm, sigmas_b = sigmas_ab, n_b = n_obs, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/I1_rmdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = TRUE,
                        rand.seed = rand.seed, included_stats = c("rmdm"))

# Plot 1A: 2D error rate of rMDM < rmu in mu space
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_rmdm_lt_rmu_dm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile() + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 2*0.05), limits=c(0, 2*0.05)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_1a rmdm 2D coverage error, mu space.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 

# Plot 1B: Heatmap of tested coverage error of rMDM < rmu in mu space
#------------------------------------------------------------------------------#
# Convert from matrix to data frame
df_err_test <- cbind(sigma = sigmas_dm, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
df_err_test$mu <- as.numeric(df_err_test$mu)
df_err_test$sigma <- as.numeric(df_err_test$sigma)
df_err_test$z <- factor(df_err_test$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df_err_test, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_1b rmdm 2D coverage error test, mu space.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)


# Plot 1C: Line plot of location of error boundaries or rmdm with mu/sigma space
#                                                                              #
source("R/coverage_error_toolbox.R")
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
mus_dm <- seq(.1, 2, by = .1)
sigmas_dm <- seq(.01, .5, by = .025)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_ab = sigmas_dm/sqrt(2/n_obs)
n_samples <- 1e3
df_crit_mu_ov_sigma <- 
  locate_bidir_binary_thresh(ind_var = "rmdm", pop_var = "rmu_dm", mus_a, sigmas_a = sigmas_ab, n_a = n_obs, 
                             mus_b = mus_a + mus_dm, sigmas_b = sigmas_ab, n_b = n_obs, mu_vsigmas_dm = NA, alphas = 0.05, n_samples, 
                             temp_path = paste(fig_path, "/I1_mdm_Error_sigma vs sigma.rds",sep=""), 
                             overwrite = FALSE, is_parallel_proc = FALSE, raw_error = FALSE, rel_error = TRUE)


# Plot of each boundary separate
df_crit_mu_ov_sigma$merge = paste(df_crit_mu_ov_sigma$er, df_crit_mu_ov_sigma$side)
gg <- ggplot(data=df_crit_mu_ov_sigma, aes(x=sigma, y=critical_mu,
                                          shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  ylab(expression(sigma[DM])) + 
  xlab(expression(mu[DM])) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
# scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  "))
gg
save_plot(paste(fig_path, "/", fig_num, "_3c rmdm error boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)






# Row 2: coverage error rate of rmdm in mu/sigma space
#
#______________________________________________________________________________#

# Coverage error simulations for mu space  
n_obs = 50
sigmas_dm <- seq(.1, 5, by = .1)
mu_vsigmas_dm <- seq(-3, 3, by = .1)

# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_obs)
sigmas_b = sigmas_a
n_samples <- 1e3

# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  quant_coverage_errors(mus_a = 100,  sigmas_a = sigmas_a, n_a = n_obs, 
                        mus_b = NA, sigmas_b = sigmas_b, n_b = n_obs, 
                        mu_vsigmas_dm = mu_vsigmas_dm, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/I1_mdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = TRUE,
                        included_stats = c("rmdm"))

# Plot 2A: error rate of rMDM < rmu in mu/sigma space
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_rmdm_lt_rmu_dm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile() + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 2*0.05), limits=c(0, 2*0.05)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_2a rmdm 2D coverage error, mu.sigma space.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 

# Plot 2B: Heatmap of tested coverage error of rMDM < rmu in mu/sigma space
#------------------------------------------------------------------------------#
# Convert from matrix to data frame
df_err_test <- cbind(sigma = sigmas_dm, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
df_err_test$mu <- as.numeric(df_err_test$mu)
df_err_test$sigma <- as.numeric(df_err_test$sigma)
df_err_test$z <- factor(df_err_test$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df_err_test, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_2b rmdm 2D coverage error test, mu.sigma space.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)







# Row 3: coverage error rate of mdm in mu/sigma space with more range and sampling
#
#______________________________________________________________________________#

# Coverage error at larger mu_dms
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
mus_dm <- seq(-100, 100, by = 2)
sigmas_dm <- seq(.01, .5, by = .025)
# mus_dm <- seq(5, 10, by = .1)
# sigmas_dm <- seq(.1, .2, by = .02)

# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_ab = sigmas_dm/sqrt(2/n_obs)
n_samples <- 1e3
mu_ov_sigmas = NULL
# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  quant_coverage_errors(mus_a = mus_a, sigmas_a = sigmas_ab, n_a = n_obs, 
                        mus_b = mus_a + mus_dm, sigmas_b = sigmas_ab, n_b = n_obs, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/I3_rmdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = TRUE,
                        rand.seed = rand.seed,
                        included_stats = c("rmdm"))

# Plot 3: 2D error rate of rMDM < rmu in mu space
#------------------------------------------------------------------------------
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_rmdm_lt_rmu_dm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile() + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 2*0.05), limits=c(0,0.13), values = c(0,0.38,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "3a rmdm coverage error extended mu_ov_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 

# Convert from matrix to data frame
df_err_test <- cbind(sigma = sigmas_dm, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
df_err_test$mu <- as.numeric(df_err_test$mu)
df_err_test$sigma <- as.numeric(df_err_test$sigma)
df_err_test$z <- factor(df_err_test$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df_err_test, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_3b rmdm coverage error extended mu_ov_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 
























































# 3rd ROW: Error rates in mu/sigma space
#       Heatmap error rate of rmdm < rmu in mu/sigma space
#       Heatmap tested coverage error of rmdm < rmu in mu/sigma space
#       Line plot of location of error boundaries or rmdm with mu/sigma space 
#_______________________________________________________________________________
# mus=NULL
# n_samples <- 1e3
# n_obs <- 50
# sigmas <- seq(.1, 5, by = .1)
# mu_ov_sigmas <- seq (-.5, .5, by=0.01)
# set.seed(rand.seed)


# First Row
# Coverage error simulations for mu space  
n_obs = 50
sigmas_dm <- seq(.1, 5, by = .1)
mu_vsigmas_dm <- seq(-3, 3, by = .1)


# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_obs)
sigmas_b = sigmas_a
n_samples <- 1e3

df_results <- 
  quant_coverage_errors(mus_a = 0.1,  sigmas_a = sigmas_a, n_a = n_obs, 
                        mus_b = NA, sigmas_b = sigmas_b, n_b = n_obs, 
                        mu_vsigmas_dm = mu_vsigmas_dm, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/rmdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE)


# # Run simulations calculating error of mdm with mu and sigma swept
# df_results <- quant_coverage_errors(NULL, sigmas, n_samples, n_obs,
#                                 paste(fig_path, "/mdm_Error_2D_mu_over_sigma_vs_sigma.rds", sep=""),
#                                 mu_ov_sigmas, overwrite = overwrite,mus_a=100, sigmas_a=1)

# Plot 9: Heatmap error rate of rmdm < rmu in mu/sigma space
#_______________________________________________________________________________
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_rmdm_lt_rmu_dm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
# grid_slopes <- slopes_by_rowcol(df_results$mean_rmdm_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM]/sigma[DM])) + ylab(expression(sigma[DM])) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 0.1), limits=c(0,0.1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_3a rmdm error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 

# Plot 10: Heatmap tested coverage error of rmdm < rmu in mu/sigma space
#_______________________________________________________________________________
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas)*length(mu_ov_sigmas)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas)*length(mu_ov_sigmas))))) %>% 
  gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM]/sigma[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_3b rmdm error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)


# Plot 11: Line plot of location of error boundaries or rmdm with mu/sigma space
#                                                                              #
#______________________________________________________________________________#
n_samples <- 1e3
n_obs <- 50
mus <- NULL
sigmas <- seq(1.5, 5, by = .1); 
mu_ov_sigmas <- seq (0.10, 0.40, by=0.001)

source("R/coverage_error_toolbox.R")

df_crit_mu_ov_sigma <- 
  locate_bidir_binary_thresh(ind_var = "rmdm", pop_var = "rmu_dm", mus=NULL, sigmas = sigmas, 
                             n_samples = n_samples, n_obs = n_obs, 
                             temp_path = paste(fig_path, "/mdm_Error_mu_over_sigma_vs_sigma.rds",sep=""), 
                             mu_ov_sigmas = mu_ov_sigmas, rand.seed = rand.seed,
                             overwrite = TRUE,  mus_a = 100, sigmas_a = 1)
# Plot of each boundary separate
df_crit_mu_ov_sigma$merge = paste(df_crit_mu_ov_sigma$er, df_crit_mu_ov_sigma$side)
gg <- ggplot(data=df_crit_mu_ov_sigma,aes(x=sigma, y=critical_mu_over_sigma,
                                          shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  ylab(expression(sigma[DM])) + 
  xlab(expression(mu[DM]/sigma[DM])) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
# scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  "))
gg
save_plot(paste(fig_path, "/", fig_num, "_3c rmdm error boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)











# Row 3: coverage error rate of mdm in mu/sigma space with more range and sampling
#______________________________________________________________________________#

# Coverage error at larger mu_dms
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
mus_dm <- seq(-50, 50, by = 2)
sigmas_dm <- seq(.01, .5, by = .025)
# mus_dm <- seq(5, 10, by = .1)
# sigmas_dm <- seq(.1, .2, by = .02)

# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_ab = sigmas_dm/sqrt(2/n_obs)
n_samples <- 1e3
mu_ov_sigmas = NULL
# Run simulations calculating error of mdm with mu and sigma swept
time1 <- Sys.time()
df_results <- 
  quant_coverage_errors(mus_a = mus_a, sigmas_a = sigmas_ab, n_a = n_obs, 
                        mus_b = mus_a + mus_dm, sigmas_b = sigmas_ab, n_b = n_obs, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "rmdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = TRUE,
                        rand.seed = rand.seed)
time2 <- Sys.time()
time2 - time1

# Plot 3: 2D error rate of rMDM < rmu in mu space
#------------------------------------------------------------------------------
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_rmdm_lt_rmu_dm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu); df$sigma <- as.numeric(df$sigma)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + geom_tile() + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) + theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       breaks = c(0, 0.05, 2*0.05), limits=c(0,0.13), values = c(0,0.38,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "3a rmdm coverage error extended mu_ov_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 

# Convert from matrix to data frame
df_err_test <- cbind(sigma = sigmas_dm, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
df_err_test$mu <- as.numeric(df_err_test$mu)
df_err_test$sigma <- as.numeric(df_err_test$sigma)
df_err_test$z <- factor(df_err_test$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df_err_test, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_3b rmdm coverage error extended mu_ov_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2, base_asp = 3, base_width = 2, dpi = 600) 





