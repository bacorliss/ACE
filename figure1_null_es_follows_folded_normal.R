

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

# 1A Show examples of histograms of normals transformed into 
set.seed(rand_seed)
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
  theme_classic(base_size=8) + theme(legend.position="none")
p_1a
# Export to TIF
save_plot("figure/figure_1a_example_dists.tiff", p_1a, ncol = 1, nrow = 1, base_height = 1.75,
          base_asp = 3, base_width = 3, dpi = 600) 




# Show that absolute of normal samples generated under null hypothesis follow folded normal
#-------------------------------------------------------------------------
# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data
set.seed(rand_seed)
n_trials <- 1e4
mus <- rep(0, n_trials)#mrunif(1e4, min = -5, max = 5)
n_obs <- 1000

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
  geom_violin(scale="width", kernel="gaussian") +
  xlab("") + ylab("P-Val.") + geom_boxplot(width=0.1,outlier.size=1) +
  theme_classic(base_size=8) + theme(legend.position="none",axis.title.x=element_blank())
p_1b
# Export to TIF
save_plot("figure/figure_1b_dist_comparison.tiff", p_1b, ncol = 1, nrow = 1, base_height = .75,
          base_asp = 3, base_width = 3, dpi = 600)


# Show % of trials that have significant p-values
# Under null disctribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_value_norm2abs_norm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_value_norm2fnorm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_value_abs_norm2fnorm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)

comp_levels <- c("N : |N|","N : F[N]","|N| : F[N]")
df_1c <- tibble(d = factor(comp_levels,levels = comp_levels, ordered = TRUE), 
       estimate = c(bin_n_an$estimate,bin_n_fn$estimate,bin_an_fn$estimate),
       lcl = c(bin_n_an$conf.int[1],bin_n_fn$conf.int[1],bin_an_fn$conf.int[1]),
       ucl = c(bin_n_an$conf.int[2],bin_n_fn$conf.int[2],bin_an_fn$conf.int[2]),
       p_value = c(bin_n_an$p.value,bin_n_fn$p.value,bin_an_fn$p.value),
       y_min = c(.997,.997,0.04),y_max = c(1, 1, 0.06),
       y_horz_line=c(1,1,0.05))


p_1c <- ggplot(df_1c, aes(x=d, y=estimate)) + 
  geom_point(size=1) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1) +
  facet_wrap(~d,scales="free", strip.position = 'bottom') +
  xlab("") +  ylab("Trials p<\u03B1 ") +
  theme_classic(base_size=8) + theme(legend.position="none",strip.text.x = element_blank(),
                                     axis.ticks.x=element_blank()) + #axis.text.x=element_blank()
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) 
p_1c  
save_plot("figure/figure_1c_dist_comparison.tiff", p_1c, ncol = 1, nrow = 1, base_height = 1,
          base_asp = 3, base_width = 2.5, dpi = 600)

# Plot hisogram, of KS values and show that they are derived from a uniform distribution
# Stats
chi_test <- chisq.test(hist(p_value_abs_norm2fnorm, breaks=20)$counts/20, simulate.p.value = TRUE)
chi_test$p.value
# Plotting
p_1c2 <- ggplot(tibble(x=p_value_abs_norm2fnorm), aes(x=x)) +
  geom_histogram(aes(y=stat(width*density)), bins=20, fill="grey") +
  theme_classic(base_size=8) + theme(legend.position="none") +
  xlab("KS P-Values, |N| : F[N]") + ylab("Probablitly") +
  geom_hline(yintercept=1/20, linetype=1, color="black",alpha=0.5, size=.5)
p_1c2             
save_plot("figure/figure_1c2_dist_comparison.tiff", p_1c2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600)
  





# Plot a single line comparing output between FN and |N|
# Plot difference between P(N) and P(F[N]) over mu
sort_x_abs_norm = sort(x_abs_norm[1,])
sort_x_fnorm = sort(x_fnorm[1,])

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
fit <- unname(lm(formula = sort_x_abs_norm ~ 0 +sort_x_fnorm)$coefficients[1])

# Define data for plotting
df_1d = tibble(x=sort_x_fnorm, y=sort_x_abs_norm)
p_1d <- ggplot(df_1d, aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab("Folded Normal Sample") + ylab("Absolute Normal Sample   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1d
save_plot("figure/figure_1d_qq_plot.tiff", p_1d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)







# Fit a line for each of the simulations
lin_model <- function(x,y) lm(formula = y ~ 0 + x)$coefficients[1];
fits <- unname(sapply(1:dim(x_fnorm)[1], function(i) lin_model(sort(x_fnorm[i,]), sort(x_abs_norm[i,])),simplify = TRUE ))


cl_fits <- quantile(fits, c(.025, .975)) 
df_fits <- tibble(x=as.factor("|N|:F[N]"),y=mean(fits))

p_1d2 <- ggplot(data=df_fits, aes(x=x,y=y)) +
  geom_point() +
  geom_hline(yintercept=1, linetype=2, color="red", alpha=.2)+
  geom_linerange(aes(ymin = cl_fits[1], ymax = cl_fits[2])) +
  ylab("Linear Regr. Coeff.") + xlab("") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1d2
save_plot("figure/figure_1d2_linear_slopes.tiff", p_1d2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1, dpi = 600)





# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
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
            y = rep(gglabel_height(df_results$sample_mean,1,.1), length(pop_mus)), label=df_mcomp$prox_mean_sig_str), 
            size=sig_char_size) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_mus),
            y = rep(gglabel_height(df_results$sample_mean,3,.1), length(pop_mus)), label=df_mcomp$prepost_mean_sig_str)) +
  xlab("Population Mean") + ylab("Sample Mean") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1e
save_plot("figure/figure_1e_mean_changes_sample_mean.tiff", p_1e, ncol = 1, nrow = 1, base_height = norm2fnorm_plot_height,
          base_asp = 3, base_width = 3, dpi = 600)

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
save_plot("figure/figure_1f_mean_changes_sample_sd.tiff", p_1f, ncol = 1, nrow = 1, base_height = norm2fnorm_plot_height,
          base_asp = 3, base_width = 3, dpi = 600)





# With fixed mean, change the sd from .1:.1:2, calculated ffolded mean and sd, and show they vary
# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
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
save_plot("figure/figure_1g_sd_changes_sample_mean.tiff", p_1g, ncol = 1, nrow = 1, base_height = norm2fnorm_plot_height,
          base_asp = 3, base_width = 3, dpi = 600)

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
save_plot("figure/figure_1h_sd_changes_sample_sd.tiff", p_1h, ncol = 1, nrow = 1, base_height = norm2fnorm_plot_height,
          base_asp = 3, base_width = 3, dpi = 600)


# With fixed mean, change the sd from .1:.1:2, calculated ffolded mean and sd, and show they vary
# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
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
  geom_text(data=data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),  
                                                         y = rep(gglabel_height(df_results$sample_mean,1,.05), length(pop_sigmas)), label=df_mcomp$prox_mean_sig_str)) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),
                                                             y = rep(gglabel_height(df_results$sample_mean,6,.05), length(pop_sigmas)), label=df_mcomp$prepost_mean_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample Mean") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1i1
save_plot("figure/figure_1i1_sd_changes_sample_mean.tiff", p_1i1, ncol = 1, nrow = 1, base_height = norm2fnorm_plot_height,
          base_asp = 3, base_width = 3, dpi = 600)

# Sample sd versus population mean
p_1i2 <- ggplot(df_results, aes(x=pop_sigma, y=sample_sd)) + 
  geom_violin(aes(fill=distr),color="black",lwd = .1, scale="width", position=position_dodge(0), alpha=0.2) + 
  geom_boxplot(aes(fill=distr), outlier.shape = NA,lwd = .1, position=position_dodge(0), alpha=0, width=.3) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),
            y = rep(gglabel_height(df_results$sample_sd,1,.05), length(pop_sigmas)), label=df_mcomp$prox_sd_sig_str)) +
  geom_text(data = data.frame(), size = sig_char_size, aes(x = factor(pop_sigmas),
            y = rep(gglabel_height(df_results$sample_sd,6,.05), length(pop_sigmas)), label=df_mcomp$prepost_sd_sig_str)) +
  coord_cartesian(clip = 'off') +
  xlab("Population SD") + ylab("Sample SD") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1i2
save_plot("figure/figure_1i2_sd_changes_sample_sd.tiff", p_1i2, ncol = 1, nrow = 1, base_height = norm2fnorm_plot_height,
          base_asp = 3, base_width = 3, dpi = 600)


