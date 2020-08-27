



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
n_samples <- 1e2

# debug(stats_param_sweep)
# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma.rds",
                                overwrite=overwrite, mus_a=100, sigmas_a=1)

# Assemble results into square matrix
grid_slopes <- rowcol_slopes(df_results$mean_diff_mmd_mu, sigmas, mus)


# Difference MMD from Mu
#------------------------------------------------------------------------------
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rmmd_error_rate)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
grid_slopes <- rowcol_slopes(df_results$mean_rmmd_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
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
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a1 rmmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 
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
  xlab(expression(sigma)) + ylab(expression(paste("Slope By Row  (ER ~ ",abs(~mu*phantom(.))," )"))) + theme_classic(base_size=8)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a3 rmmd error rate 2D_vs_sigmas.tiff",sep=""),
          gg, base_height = 2.2, base_asp = 3, base_width = 2.2, dpi = 600) 




# 2D visualization of hypothesized coverage error of MMD in mu space
#                                                                              #
#______________________________________________________________________________#
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$p_val_mmd_eq_zero > 0.05/(length(sigmas)*length(mus)),
  df_results$p_val_mmd_eq_alpha > 0.05/(length(sigmas)*length(mus))))) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) +
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






# 2D visualization of hypothesized coverage error of MMD in mu/sigma space
#                                                                              #
#______________________________________________________________________________#
sigmas <- seq(.1, 5, by = .1)
mu_ov_sigmas <- seq (-.5, .5, by=0.01)
n_obs <- 50
set.seed(rand.seed)
# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(NULL, sigmas, n_samples, n_obs,
                           "temp/mmd_Error_2D_mu_over_sigma_vs_sigma.rds", mu_ov_sigmas)
# Difference MMD from Mu
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$p_val_mmd_eq_zero > 0.05/(length(sigmas)*length(mu_ov_sigmas)),
  df_results$p_val_mmd_eq_alpha > 0.05/(length(sigmas)*length(mu_ov_sigmas))))) %>% 
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
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_b5 mmd error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)






# Identify location of coverage error boundaries with mu space
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
sigmas <- seq(.1, 5, by = .1); 
mus <- seq (.1/5, 2, by=0.02)
mu_ov_sigmas = NULL

df_crit_mu <- locate_bidir_binary_thresh(mus = mus, sigmas = sigmas, n_obs = n_obs, 
                                         temp_name = "mmd_Error_mu_vs_sigma.rds",
                                         mu_ov_sigmas = mu_ov_sigmas, rand.seed = rand.seed)
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
  xlab(expression(sigma)) + 
  ylab(expression(mu)) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))+
  # scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  ")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_d mmd boundaries over mu.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)

# Box and Whiskers of Coverage Error Transition Region
p <- ggplot(df_crit_mu, aes(x=er, color = side,#group = interaction(er, side), 
                            y = abs(critical_mu))) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9), outlier.shape = NA) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") + 
  ylab(expression(abs(~mu~phantom(.))))+
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
p
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_c mmd transition over mu.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(abs(critical_mu) ~ er + side, data = df_crit_mu)
summary_mu <- summary(res.aov2)
capture.output(summary_mu, file = paste("figure/F", fig_num, "/F", fig_num,
                                        "_c mmd transition over mu.txt", sep=""))


# Identify location of coverage error boundaries with mu/sigma space
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
mus <- NULL
sigmas <- seq(.1, 5, by = .1); 
mu_ov_sigmas <- seq (0.10, 0.40, by=0.001)

df_crit_mu_ov_sigma <- locate_bidir_binary_thresh(mus=NULL, sigmas = sigmas, n_obs = n_obs, 
                                      temp_name = "mmd_Error_mu_over_sigma_vs_sigma.rds", 
                                      mu_ov_sigmas = mu_ov_sigmas, rand.seed = rand.seed)
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

# Box and Whiskers of Coverage Error Transition Region
p <- ggplot(df_crit_mu_ov_sigma, aes(x=er, color = side,#group = interaction(er, side), 
                         y = abs(critical_mu_over_sigma))) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9), outlier.shape = NA) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") +  
  ylab(expression(abs(~mu*phantom(.))*phantom(.)/sigma))+
scale_color_manual(values = c("#66c2a5", "#fc8d62"))
p
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_c mmd transition over mu sigma.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(abs(critical_mu_over_sigma) ~ er + side , data = df_crit_mu_ov_sigma)
summary_mu_ov_sigma <- summary(res.aov2)
capture.output(summary_mu_ov_sigma, file = paste("figure/F", fig_num, "/F", fig_num,
                                        "_c mmd transition over mu sigma.txt", sep=""))

