

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
source("R/row_stats_toolbox.R")
# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
# Script Parameters
fig_num = "5"
fig_path = file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE, recursive=TRUE)
rand.seed <- 0
overwrite <- TRUE
yaxis_font_size <- 8



# Investigation 1: Initial coverage error of rmdm in mu space
#
#_______________________________________________________________________________
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
mus_dm <- seq(-2.5, 2.5, by = .05)#by = .025)
sigmas_dm <- seq(.01, 1, by = .01)#by = .005)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_ab = sigmas_dm/sqrt(2/n_obs)
n_samples <- 1e3
mu_ov_sigmas = NULL
# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  quant_coverage_errors(
    mus_a = mus_a, sigmas_a = sigmas_ab, n_a = n_obs, 
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
                       breaks = c(0, 0.05, .1), limits=c(0, .1)) +
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
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 
    0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 
    0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
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


# Plot 1C: Line plot of location of error boundaries or rmdm with mu space
#------------------------------------------------------------------------------#
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
mus_dm <- seq(.1, 1.2, by = .1)
sigmas_dm <- seq(.01, .5, by = .02)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_ab = sigmas_dm/sqrt(2/n_obs)
n_samples <- 1e3
df_crit_mu <- 
  locate_bidir_binary_thresh(
    ind_var = "rmdm", pop_var = "rmu_dm", mus_a, sigmas_a = sigmas_ab, n_a = n_obs, 
    mus_b = mus_a + mus_dm, sigmas_b = sigmas_ab, n_b = n_obs, 
    mu_vsigmas_dm = NA, alphas = 0.05, n_samples, 
    temp_path = paste(fig_path, "/I1_rmdm_boundary.rds",sep=""), 
    overwrite = overwrite, is_parallel_proc = TRUE, raw_error = FALSE, rel_error = TRUE)
df_slopes1 <- df_crit_mu %>% group_by(er,side) %>% summarize(pearson = cor.test(
  critical_mu, sigma,method = "pearson")$p.value)
df_slopes1$adj_pearson <- p.adjust(df_slopes1$pearson,"bonferroni")
# Plot of each boundary separate
gg <- ggplot(data=df_crit_mu, aes(x=sigma, y=critical_mu,
                                          shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  xlab(expression(sigma[DM])) + 
  ylab(expression(mu[DM])) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
# scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  "))
gg
save_plot(paste(fig_path, "/", fig_num, "_1c rmdm error boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)






# Row 2: coverage error rate of rmdm in mu/sigma space
#
#______________________________________________________________________________#

n_obs = 50
sigmas_dm <- seq(.1, 5, by = .05)
mu_vsigmas_dm <- seq(-3, 3, by = .05)

# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_obs)
sigmas_b = sigmas_a
n_samples <- 1e3

# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  quant_coverage_errors(mus_a = 100,  sigmas_a = sigmas_a, n_a = n_obs, 
                        mus_b = NA, sigmas_b = sigmas_b, n_b = n_obs, 
                        mu_vsigmas_dm = mu_vsigmas_dm, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/I2_rmdm_Error_2D_mu_vs_sigma.rds",sep=""),
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
  xlab(expression(mu[DM]/sigma[DM])) + ylab(expression(sigma[DM])) + theme_classic(base_size=8) +
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
  df_results$pval_err_eq_zero_abs_rmdm_lt_rmu_dm > 
    0.05/(length(sigmas_dm)*length(mus_dm)),
  df_results$pval_err_eq_alpha_abs_rmdm_lt_rmu_dm > 
    0.05/(length(sigmas_dm)*length(mus_dm))))) %>% gather(mu, z, -sigma)
df_err_test$mu <- as.numeric(df_err_test$mu)
df_err_test$sigma <- as.numeric(df_err_test$sigma)
df_err_test$z <- factor(df_err_test$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df_err_test, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM]/sigma[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "red", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste(fig_path, "/", fig_num, "_2b rmdm 2D coverage error test, mu.sigma space.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)

# Plot 2C: Line plot of location of error boundaries or rmdm with mu/sigma space
#------------------------------------------------------------------------------#
# Coverage error simulations for mu space  
n_obs = 100
mus_a = 10
# mus_dm <- seq(.1, 1.2, by = .1)
sigmas_dm <- seq(.11, .5, by = .01)
mu_vsigmas_dm <- seq(1.3, 2, by = .025)
# Spread sigma_dm across sigma_a and sigma_b equally
sigmas_a = sigmas_dm/sqrt(2/n_obs)
sigmas_b = sigmas_a
n_samples <- 1e3

df_crit_mu_ov_sigma <- 
  locate_bidir_binary_thresh(
    ind_var = "rmdm", pop_var = "rmu_dm", mus_a, sigmas_a = sigmas_a, n_a = n_obs,
    mus_b = NA, sigmas_b = sigmas_b, n_b = n_obs, mu_vsigmas_dm = mu_vsigmas_dm, alphas = 0.05, n_samples,
    temp_path = paste(fig_path, "/I2_rmdm_error_boundary_mu-sigma_vs_sigma.rds",sep=""), 
    overwrite = overwrite, is_parallel_proc = TRUE, raw_error = FALSE, rel_error = TRUE)
df_crit_mu_ov_sigma$merge = paste(df_crit_mu_ov_sigma$er, df_crit_mu_ov_sigma$side)
# Test for relationship between border and mu/sigma
df_slopes2 <- df_crit_mu_ov_sigma %>% group_by(er,side) %>% summarize(pearson = cor.test(
  critical_mu_over_sigma, sigma,method = "pearson")$p.value, x_bar = mean(critical_mu_over_sigma))
# df_slopes2$adj_pearson <- p.adjust(df_slopes2$pearson,"bonferroni")
# df_slopes2 <- df_crit_mu_ov_sigma %>% group_by(er,side) %>% 
#   summarize(ci = confint(lm(critical_mu_over_sigma~sigma),'sigma',level=0.95))
# Plot of each boundary separate
gg <- ggplot(data=df_crit_mu_ov_sigma, aes(x=sigma, y=critical_mu_over_sigma,
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
save_plot(paste(fig_path, "/", fig_num, "_2c rmdm error boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)

