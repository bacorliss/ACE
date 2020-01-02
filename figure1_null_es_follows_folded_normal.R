

# The magntiude of effect sizefrom normal data follows the folded normal 
# distribution and converges to the normal distribution as mu>>sigma




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
set.seed(0)
dist_text_size=3;
equ_font_size=8;

# Initial parameters
nsamples = 1e6
ngroups = 2
set.seed(0)
dist_text_size=3;
equ_font_size=8;





# 1A Show examples of histograms of normals transformed into 
sample_1a = rnorm(nsamples, mean = 0, sd = 1)
df_1a = rbind(tibble(group = 2, mu = "mu==0", x = sample_1a    ), tibble(group = 1, mu = "mu==0", x = abs(sample_1a)    ),
              tibble(group = 2, mu = "mu==1", x = sample_1a + 1), tibble(group = 1, mu = "mu==1", x = abs(sample_1a + 1)),
              tibble(group = 2, mu = "mu==1.5", x = sample_1a + 1.5), tibble(group = 1, mu = "mu==1.5", x = abs(sample_1a + 1.5)),
              tibble(group = 2, mu = "mu==3", x = sample_1a + 3), tibble(group = 1, mu = "mu==3", x = abs(sample_1a + 3))
              )
df_1a$group <- as.factor(df_1a$group)
df_1a$mu    <- as.factor(df_1a$mu)

df_1a %>%
  group_by(group,mu) %>%
  summarise(mean = mean(x), n = n())

p_1a <- ggplot(df_1a, aes(x = x)) +
  geom_histogram(aes(fill = group,y=..density..), binwidth = .2, boundary = 0, position = "identity", alpha = 0.5) + 
  facet_grid(. ~ mu, scales = "free_x", labeller=label_parsed) + ylab("Probability") +
  theme_classic(base_size=8)
p_1a
# Export to TIF
save_plot("figure/figure_1a_example_dists.tiff", p_1a, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.5, dpi = 600) 




# Show that absolute of normal samples generated under null hypothesis follow folded normal
#-------------------------------------------------------------------------

# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data                |
n_trials <- 1e3
mus <- rep(0, n_trials)#mrunif(1e4, min = -5, max = 5)
n_obs <- 50

# Generate random samples based on random mu values
x_norm <- t(sapply(mus, function (x) rnorm(n_obs, mean = x, sd = 1),simplify = TRUE))
x_abs_norm <- t(sapply(mus, function (x) abs(rnorm(n_obs, mean = x, sd = 1)),simplify = TRUE))
x_fnorm <- t(sapply(mus, function (x) rfoldnorm(n_obs, mean = x, sd = 1),simplify = TRUE))

# Test if norm is different from abs_norm
p_value_norm2abs_norm <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_abs_norm[i,]), alternative = "two.sided", exact = TRUE, 
                                tol=1e-8, simulate.p.value=FALSE)$p)

# Test if norm is different from folded norm
p_value_norm2fnorm <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                             as.vector(x_fnorm[i,]), alternative = "two.sided", exact = TRUE, 
                             tol=1e-8, simulate.p.value=FALSE)$p)

# Test if abs_norm is different from folded norm
p_value_abs_norm2fnorm <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_abs_norm[i,]),
                                 as.vector(x_fnorm[i,]), alternative = "two.sided", exact = TRUE, 
                                 tol=1e-8, simulate.p.value=FALSE)$p)
df_1b <- rbind(tibble(d = factor("N : |N|"), x = p_value_norm2abs_norm),
               tibble(d = factor("N : F[N]"), x = p_value_norm2fnorm), 
               tibble(d = factor("|N| : F[N]"), x = p_value_abs_norm2fnorm))

p_1b <- ggplot(df_1b, aes(x=d, y=x)) + 
  geom_violin(scale="width",bw=.05,kernel="gaussian") +
  xlab("Sample Comparisons") + ylab("K.S. P Value") + geom_boxplot(width=0.1,outlier.shape = NA) +
  theme_classic(base_size=8)
p_1b
# Export to TIF
save_plot("figure/figure_1b_dist_comparison.tiff", p_1b, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600)




# Show % of trials that have significant p-values
# Under null disctribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_value_norm2abs_norm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_value_norm2fnorm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_value_abs_norm2fnorm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)

df_1c <- tibble(d = factor(c("N : |N|","N : F[N]","|N| : F[N]")), 
       estimate = c(bin_n_an$estimate,bin_n_fn$estimate,bin_an_fn$estimate),
       lcl = c(bin_n_an$conf.int[1],bin_n_fn$conf.int[1],bin_an_fn$conf.int[1]),
       ucl = c(bin_n_an$conf.int[2],bin_n_fn$conf.int[2],bin_an_fn$conf.int[2]),
       p_value = c(bin_n_an$p.value,bin_n_fn$p.value,bin_an_fn$p.value))

p_1c <- ggplot(df_1c, aes(x=d, y=estimate)) + 
  geom_hline(yintercept = 0.05, linetype="dashed", color="red", size=.5,alpha=.2) +
  geom_hline(yintercept = 1, linetype="dashed", color="black", size=.5,alpha=.2) +
  geom_point(size=.1) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2) +
  xlab("Sample Comparisons") +  ylab("Fraction Trials p<a  ") +
  theme_classic(base_size=8)
p_1c  
save_plot("figure/figure_1c_dist_comparison.tiff", p_1c, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)




# Plot difference between P(N) and P(F[N]) over mu
sort_x_abs_norm = sort(x_abs_norm[1,])
sort_x_fnorm = sort(x_fnorm[1,])

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
fit <- lm(formula = sort_x_abs_norm ~ 0 +sort_x_fnorm)
confint(fit,level=0.95)
# Define data for plotting
df_1d = tibble(x=sort_x_fnorm, y=sort_x_abs_norm)

p_1d <- ggplot(df_1d, aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab("Folded Normal Sample") + ylab("Absolute Normal Sample   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_classic(base_size=8)
p_1d
save_plot("figure/figure_1d_qq_plot.tiff", p_1d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)






# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
pop_mus =  seq(-3,3, by = 1)
pop_sigmas = 1
n_obs = 100
n_samples = 1e4
# Run Simulation, calculate mean and sd for each sample
df_results <- rfnorm_swept_param(pop_mus,pop_sigmas,n_samples, n_obs)

# Statistics of groups
# Compute the analysis of variance


# summary(aov(sample_mean ~ pop_mu, data = df_results))
# pairwise.t.test(df_results$sample_mean, df_results$pop_mu, p.adj = "bonf")
# pairwise.t.test(df_results$sample_sd, df_results$pop_mu, p.adj = "bonf")
# 
# prox_sd_sig <- rep('#', length(pop_mus))
# prox_mean_sig <- rep('#', length(pop_mus))
# prepost_sd_sig <- rep('\u2020', length(pop_mus))
# prepost_mean_sig <- rep('\u2020', length(pop_mus))
# df_mcomp = tibble(pop_mu=pop_mus, prox_sd_sig=prox_sd_sig, prox_mean_sig=prox_mean_sig,
#                   prepost_sd_sig,prepost_mean_sig)

# Calculate adjancent and pairwise stats of results
df_mcomp <- norm_fnorm_stats(df_results, "pop_mu") 

# Sample mean versus population mean
p_1e <- ggplot(df_results, aes(x=pop_mu, y=sample_mean)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .2, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .2, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), aes(x = factor(pop_mus),  
            y = rep(gglabel_height(df_results$sample_mean,1,.1), length(pop_mus)), label=df_mcomp$prox_mean_sig_str), 
            size=2) +
  geom_text(data = data.frame(), size = 2, aes(x = factor(pop_mus),
            y = rep(gglabel_height(df_results$sample_mean,3,.1), length(pop_mus)), label=df_mcomp$prepost_mean_sig_str)) +
  xlab("Population Mean") + ylab("Sample Mean") +
  theme_classic(base_size=8)
p_1e
save_plot("figure/figure_1e_mean_changes_sample_mean.tiff", p_1e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600)

# Sample sd versus population mean
p_1f <- ggplot(df_results, aes(x=pop_mu, y=sample_sd)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .2, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .2, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), aes(x = factor(pop_mus),  
            y = rep(gglabel_height(df_results$sample_sd,1,.1), length(pop_mus)), label=df_mcomp$prox_sd_sig_str), 
            size=2) +
  geom_text(data = data.frame(), size = 2, aes(x = factor(pop_mus),
            y = rep(gglabel_height(df_results$sample_sd,3,.1), length(pop_mus)), label=df_mcomp$prepost_sd_sig_str)) +  
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample SD") +
  theme_classic(base_size=8)
p_1f
save_plot("figure/figure_1f_mean_changes_sample_sd.tiff", p_1f, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600)





# With fixed mean, change the sd from .1:.1:2, calculated ffolded mean and sd, and show they vary
# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
pop_mus = 0
pop_sigmas = seq(.5,5, by = .5)
n_obs = 100
n_samples = 1e4
# Run Simulation, calculate mean and sd for each sample
df_results <- rfnorm_swept_param(pop_mus,pop_sigmas,n_samples = 1e4, n_obs = 100)
  
# # Statistics of groups
# # Compute the analysis of variance
# summary(aov(sample_mean ~ pop_sigma, data = df_results))
# summary(aov(sample_sd ~ pop_sigma, data = df_results))
# pw_sample_mean <- pairwise.t.test(df_results$sample_mean, df_results$pop_sigma, p.adj = "bonf")
# pw_sample_sd <- pairwise.t.test(df_results$sample_sd, df_results$pop_sigma, p.adj = "bonf")
# 
# prox_sd_sig <- rep('#', length(pop_sigmas))
# prox_mean_sig <- rep('#', length(pop_sigmas))
# prepost_sd_sig <- rep('\u2020', length(pop_sigmas))
# prepost_mean_sig <- rep('\u2020', length(pop_sigmas))
# df_mcomp = tibble(pop_sigma=pop_sigmas, prox_sd_sig=prox_sd_sig, prox_mean_sig=prox_mean_sig,
#                   prepost_sd_sig,prepost_mean_sig)

# Calculate adjancent and pairwise stats of results
df_mcomp <- norm_fnorm_stats(df_results, "pop_sigma") 


# Sample mean versus population mean
p_1g <- ggplot(df_results, aes(x=pop_sigma, y=sample_mean))+ 
  geom_violin(aes(fill=distr),color="black",lwd = .2, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .2, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data=data.frame(), aes(x = factor(pop_sigmas),  
            y = rep(gglabel_height(df_results$sample_mean,1,.05), length(pop_sigmas)), label=df_mcomp$prox_mean_sig_str), 
            size=2) +
  geom_text(data = data.frame(), size = 2, aes(x = factor(pop_sigmas),
            y = rep(gglabel_height(df_results$sample_mean,6,.05), length(pop_sigmas)), label=df_mcomp$prepost_mean_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample Mean") +
  theme_classic(base_size=8)
p_1g
save_plot("figure/figure_1g_sd_changes_sample_mean.tiff", p_1g, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600)

# Sample sd versus population mean
p_1h <- ggplot(df_results, aes(x=pop_sigma, y=sample_sd)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .2, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .2, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data = data.frame(), size = 2, aes(x = factor(pop_sigmas),
            y = rep(gglabel_height(df_results$sample_sd,1,.05), length(pop_sigmas)), label=df_mcomp$prox_sd_sig_str)) +
  geom_text(data = data.frame(), size = 2, aes(x = factor(pop_sigmas),
             y = rep(gglabel_height(df_results$sample_sd,6,.05), length(pop_sigmas)), label=df_mcomp$prepost_sd_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample SD") +
  theme_classic(base_size=8)
p_1h
save_plot("figure/figure_1h_sd_changes_sample_sd.tiff", p_1h, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600)




