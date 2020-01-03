


# The absolute value transforms alters normal data in a nonlinear maner, and the population parameters mu[f] and
# sigma f of a folded normal are dependent on each other.


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

source("R/stat_helper.r")

# choose colors for plotting
color_pal = brewer.pal(4, "Set1")


# Parameters
rand_seed <- 0
dist_text_size=3;
equ_font_size=8;



# Initial parameters
nsamples = 1e4
ngroups = 2
set.seed(0)
dist_text_size=3
equ_font_size=8
sig_char_size = 3
norm2fnorm_plot_height = 1.75
plot_height = 1.25

fig_basename="sf_2"




# Show how folded sample mean is changed as mu increases from zero --------------------------------------
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
set.seed(rand_seed)
pop_mus =  seq(-4.5,4.5, by = 1.5)
pop_sigmas = 1
n_obs = 100
n_samples = 1e3
# Run Simulation, calculate mean and sd for each sample
df_results <- rfnorm_swept_param(pop_mus,pop_sigmas,n_samples, n_obs)

# Calculate adjancent and pairwise stats of results
df_mcomp <- norm_fnorm_stats(df_results, "pop_mu") 

# Sample mean versus population mean
p_1e <- ggplot(df_results, aes(x=pop_mu, y=sample_mean)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), aes(x = factor(pop_mus),  
                                   y = rep(gglabel_height(df_results$sample_mean,1,.1), length(pop_mus)), 
                                   label=df_mcomp$prox_mean_sig_str),  size=sig_char_size) +
  geom_text(data = data.frame(), size = sig_char_size, 
            aes(x = factor(pop_mus), y = rep(gglabel_height(df_results$sample_mean,3,.1), length(pop_mus)),
                label=df_mcomp$prepost_mean_sig_str)) +
  xlab("Population Mean") + ylab("Sample Mean") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1e
save_plot(paste("figure/", fig_basename, "a_mean_changes_sample_mean.tiff"), p_1e, ncol = 1, nrow = 1, 
          base_height = norm2fnorm_plot_height, base_asp = 3, base_width = 3, dpi = 600)

# Sample sd versus population mean
p_1f <- ggplot(df_results, aes(x=pop_mu, y=sample_sd)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), aes(x = factor(pop_mus),  
                                   y = rep(gglabel_height(df_results$sample_sd,1,.1), length(pop_mus)), label=df_mcomp$prox_sd_sig_str), 
            size=sig_char_size) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_mus),
                                                           y = rep(gglabel_height(df_results$sample_sd,3,.1), length(pop_mus)), label=df_mcomp$prepost_sd_sig_str)) +  
  coord_cartesian(clip = 'off') +
  xlab("Population Mean") + ylab("Sample SD") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1f
save_plot(paste("figure/",fig_basename, "b_mean_changes_sample_sd.tiff", sep = ""), p_1f, ncol = 1, nrow = 1,
          base_height = norm2fnorm_plot_height, base_asp = 3, base_width = 3, dpi = 600)




# With fixed mean, change the sd from .1:.1:2, calculated ffolded mean and sd, and show they vary
# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
set.seed(rand_seed)
pop_mus = 0
pop_sigmas = seq(.5,5.5, by = 1)
n_obs = 100
n_samples = 1e3
# Run Simulation, calculate mean and sd for each sample
df_results <- rfnorm_swept_param(pop_mus,pop_sigmas,n_samples, n_obs)
# Calculate adjancent and pairwise stats of results
df_mcomp <- norm_fnorm_stats(df_results, "pop_sigma") 

# Sample mean versus population mean
p_1g <- ggplot(df_results, aes(x=pop_sigma, y=sample_mean))+ 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),  
                                                         y = rep(gglabel_height(df_results$sample_mean,1,.05), length(pop_sigmas)), label=df_mcomp$prox_mean_sig_str)) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),
                                                           y = rep(gglabel_height(df_results$sample_mean,6,.05), length(pop_sigmas)), label=df_mcomp$prepost_mean_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample Mean") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1g
save_plot(paste("figure/", fig_basename, "c_sd_changes_sample_mean.tiff", sep = ""), p_1g, ncol = 1, nrow = 1, 
          base_height = norm2fnorm_plot_height, base_asp = 3, base_width = 3, dpi = 600)


# Sample sd versus population mean
p_1h <- ggplot(df_results, aes(x=pop_sigma, y=sample_sd)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),
                                                           y = rep(gglabel_height(df_results$sample_sd,1,.05), length(pop_sigmas)), label=df_mcomp$prox_sd_sig_str)) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),
                                                           y = rep(gglabel_height(df_results$sample_sd,6,.05), length(pop_sigmas)), label=df_mcomp$prepost_sd_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample SD") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1h
save_plot(paste("figure/", fig_basename, "d_sd_changes_sample_sd.tiff", sep = ""), p_1h, ncol = 1, nrow = 1, 
          base_height = norm2fnorm_plot_height, base_asp = 3, base_width = 3, dpi = 600)



# At large mu, folding does not change data from normal -------------------------------------------
# With fixed mean, change the sd from .1:.1:2, calculated ffolded mean and sd, and show they vary
# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
set.seed(rand_seed)
pop_mus = 10
pop_sigmas = seq(.5,5.5, by = 1)
n_obs = 100
n_samples = 1e3
# Run Simulation, calculate mean and sd for each sample
df_results <- rfnorm_swept_param(pop_mus,pop_sigmas,n_samples, n_obs)
# Calculate adjancent and pairwise stats of results
df_mcomp <- norm_fnorm_stats(df_results, "pop_sigma") 

# Sample mean versus population mean
p_1i1 <- ggplot(df_results, aes(x=pop_sigma, y=sample_mean))+ 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), size = sig_char_size, 
            aes(x = factor(pop_sigmas), y = rep(gglabel_height(df_results$sample_mean,1,.05), 
                                                length(pop_sigmas)), label=df_mcomp$prox_mean_sig_str)) +
  geom_text(data = data.frame(), size = sig_char_size, 
            aes(x = factor(pop_sigmas), y = rep(gglabel_height(df_results$sample_mean,6,.05), 
                                                length(pop_sigmas)), label=df_mcomp$prepost_mean_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample Mean") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1i1
save_plot(paste("figure/", fig_basename, "e_sd_changes_sample_mean.tiff", sep = ""), p_1i1, ncol = 1, nrow = 1,
          base_height = norm2fnorm_plot_height, base_asp = 3, base_width = 3, dpi = 600)



# Sample sd versus population mean
p_1i2 <- ggplot(df_results, aes(x=pop_sigma, y=sample_sd)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data = data.frame(), size = sig_char_size, 
            aes(x = factor(pop_sigmas), y = rep(gglabel_height(df_results$sample_sd,1,.05), 
                                                length(pop_sigmas)), label=df_mcomp$prox_sd_sig_str)) +
  geom_text(data = data.frame(), size = sig_char_size, 
            aes(x = factor(pop_sigmas), y = rep(gglabel_height(df_results$sample_sd,6,.05), 
                                                length(pop_sigmas)), label=df_mcomp$prepost_sd_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample SD") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1i2
save_plot(paste("figure/", fig_basename, "f_sd_changes_sample_sd.tiff", sep = ""), p_1i2, ncol = 1, nrow = 1,
          base_height = norm2fnorm_plot_height, base_asp = 3, base_width = 3, dpi = 600)



