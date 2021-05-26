

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
# User defined functions
source("R/mdm.R")
source("R/coverage_error_toolbox.R")


# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_z"
# Script Parameters
fig_num = "4"
fig_path = file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE,recursive = TRUE)
n_samples <- 1e3
rand.seed <- 0
overwrite <- TRUE
is_parallel_proc <- TRUE



# 2D visualization of mdm difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#

# First Row
# Coverage error simulations for mu space  
n_obs = 50
mus_dm <- seq(-2.5, 2.5, by = .1)
sigmas_dm <- seq(.01, 1, by = .02)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_obs)
sigmas_b = sigmas_a
n_samples <- 1e3
mu_vsigmas_dm = NULL


# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  quant_coverage_errors(mus_a = 100, sigmas_a = sigmas_a, n_a = n_obs, 
                        mus_b = 100 + mus_dm, sigmas_b = sigmas_b, n_b = n_obs, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/mdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE)


# Error rate of MDM < mu
#------------------------------------------------------------------------------
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_mdm_lt_mu_dm)) %>% gather(mu, z, -sigma)
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
                       limits=c(0,.1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_a3 mdm error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 


# 2D visualization of hypothesized coverage error of MDM in mu space
#                                                                              #
#______________________________________________________________________________#
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas_dm, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_mdm_lt_mu_dm > 0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_mdm_lt_mu_dm > 0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_b4 mdm error test.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)






# Row 2: error rate of mdm in mu/sigma space
#                                                                              #
#______________________________________________________________________________#


# First Row
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
                        n_samples = n_samples, out_path = paste(fig_path, "/mdm_Error_2D_mu_vs_sigma.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE)


# Error rate of MDM < mu in mu/sigma space
#------------------------------------------------------------------------------
# Convert from matrix to dataframe
df <- cbind(sigma = sigmas_dm, as_tibble(df_results$mean_err_abs_mdm_lt_mu_dm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
# grid_slopes <- slopes_by_rowcol(df_results$mean_mdm_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM]/sigma[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(0,.1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "\\", fig_num, "_3a mdm error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 

# Tested error rate MDM < mu in mu/sigma space
#-------------------------------------------------------------------------------
df <- cbind(sigma = sigmas_dm, as_tibble(error_test_codes(
  df_results$pval_err_eq_zero_abs_mdm_lt_mu_dm > 0.05/(length(sigmas_dm)*length(mu_vsigmas_dm)),
  df_results$pval_err_eq_alpha_abs_mdm_lt_mu_dm > 0.05/(length(sigmas_dm)*length(mu_vsigmas_dm))))) %>% 
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
save_plot(paste(fig_path, "\\", fig_num, "_3b mdm error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)






# Identify location of coverage error boundaries with mu space
#                                                                              #
#______________________________________________________________________________#

# First Row
# Coverage error simulations for mu space  
n_a = 50
n_b = 50
mus_dm <- seq (.02, 2, by=0.02)
sigmas_dm <- seq(.1, 1, by = .025)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_a)
sigmas_b = sigmas_a
n_samples <- 1e3


df_crit_mu <-
  locate_bidir_binary_thresh(mus_a = 100, sigmas_a = sigmas_a, n_a = n_a, 
                             mus_b = 100 + mus_dm, sigmas_b = sigmas_b, n_b = n_b,
                             mu_vsigmas_dm = NA, alphas = 0.05,
                             n_samples = n_samples, temp_path = paste(fig_path, "/mdm_Error_2D_mu_vs_sigma.rds",sep=""),
                             overwrite = overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE)
df_crit_mu$merge = paste(df_crit_mu$er, df_crit_mu$side)

# Plot of each boundary separate
df_slopes <- df_crit_mu %>% group_by(er,side) %>% summarize(pearson = cor.test(
  critical_mu, sigma,method = "pearson")$p.value)
# Plot of each boundary separate
# df_pearson$adj_pearson <- p.adjust(df_pearson$pearson,"bonferroni")
gg <- ggplot(data=df_crit_mu,aes(x=sigma, y=critical_mu,
                                 shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  xlab(expression(sigma[DM])) + 
  ylab(expression(mu[DM])) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))+
  # scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  ")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
gg
save_plot(paste(fig_path, "\\", fig_num, "_d mdm boundaries over mu.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)

# 
res.aov2 <- aov(abs(critical_mu) ~ er + side, data = df_crit_mu)
summary_mu <- summary(res.aov2)
capture.output(summary_mu, file = paste(fig_path, "\\", fig_num,
                                        "_c mdm transition over mu.txt", sep=""))




# Identify location of coverage error boundaries with mu/sigma space
#                                                                              #
#______________________________________________________________________________#


# First Row
# Coverage error simulations for mu space  
n_a = 50
n_b = 50
mus_dm <- seq (.02, 2, by=0.02)
mu_vsigmas_dm <- seq (1.5, 2.5, by=0.05)
sigmas_dm <- seq(.1, 1, by = .025)

# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_a)
sigmas_b = sigmas_a
n_samples <- 1e3

df_crit_mu_ov_sigma <- 
  locate_bidir_binary_thresh(mus_a = 100, sigmas_a = sigmas_a, n_a = n_a, 
                             mus_b = NA, sigmas_b = sigmas_b, n_b = n_b,
                             mu_vsigmas_dm = mu_vsigmas_dm, alphas = 0.05,
                             n_samples = n_samples, temp_path = paste(fig_path, "/mdm_Error_2D_mu_vs_sigma.rds",sep=""),
                             overwrite = overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE)
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
  xlab(expression(sigma[DM])) + 
  ylab(expression(mu[DM]/sigma[DM])) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
  # scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  "))
gg
save_plot(paste(fig_path, "\\", fig_num, "_d mdm boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)

res.aov2 <- aov(abs(critical_mu_over_sigma) ~ er + side , data = df_crit_mu_ov_sigma)
summary_mu_ov_sigma <- summary(res.aov2)
capture.output(summary_mu_ov_sigma, file = paste(fig_path, "\\", fig_num,
                                        "_c mdm transition over mu sigma.txt", sep=""))


