



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

source("R/mhd.r")

fig_basename = "f_3"


# Visualize T1E (mu<MHD) under the null and alt hypothesis, mu swept (sigma=1)
#
#______________________________________________________________________________#

mus <- c(seq(-10, -2, by = 1), seq(-1, 1, by = 0.05), seq(2, 10, by = 1))
sigmas <- c(1)
n_samples <- 1e4
n_obs <- 50
# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(mus), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(mus), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(mus), ncol = 3),
                mhd_true = rep(0,length(mus)))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(mus), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- sapply(rep(mus[n],n_samples), function (x)
    rnorm(n_obs, mean = mus[n], sd = 1), simplify = TRUE)
  # Calculat ethe MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[n]),sigmas[1]/sqrt(n_obs)))
  
  # Calculate the quantile of the MHD for each sample compared to pooled sample
  mhd_quantile_sample = 
     sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mean(x1)), sd(x1)/ sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(mu = mus[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
           mhd_quantile_pop = mhd_quantile_pop, 
           mhd_quantile_sample = mhd_quantile_sample)
  
  # Calculate summary stats
  df_sum$ci_95_sample[n,] = quantile(mhd_quantile_sample, probs = c(0.025,0.975))
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  df_sum$mhd_true[n] = sum(mhd>abs(mus[n]))
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[n])), n_samples, 
                                       conf.level=0.95,methods="exact"))[c(4,5,6)]
  
}
df <- do.call("rbind", list_df)

# Violin of x_bar
p_a1 <- ggplot(subset(df, mu == c( seq(-3, 3, by = 1.5))), aes(x = mu, y = x_bar)) + 
  geom_violin(scale="width", kernel="gaussian", fill="grey96", 
              color="red", draw_quantiles = c(0.95)) +
  geom_violin(scale="width", kernel="gaussian", alpha =0, 
              color="black") +
  facet_wrap(. ~ mu, scales = "free", labeller=label_parsed,ncol = length(mus)) +
  xlab(expression(mu)) + ylab(expression(bar(x))) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank())
p_a1
save_plot(paste("figure/", fig_basename, "d integration_cdf_central.tiff", 
                sep = ""), p_a1, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

# 0.95 quantile for each distribution
# MError of the 95% quantile 
p_a2 <- ggplot(df_sum, aes(x = mu, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(mu)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_a2
save_plot(paste("figure/", fig_basename, "a2 t1e_MHD_over_mu.tiff", 
                sep = ""), p_b, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# 0.95 quantile for each distribution
# MError of the 95% quantile 
p_a3 <- ggplot(subset(df_sum,mu>=-0.75 & mu<=0.75), aes(x = mu, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha=.2) +
  xlab(expression(mu)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_a3
save_plot(paste("figure/", fig_basename, "a3 t1e_MHD_over_mu_zoomed.tiff", 
                sep = ""), p_c, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)





# Visualize T1E (mu<MHD) under the null hypothesis (mu=0), sigma swept
#
#______________________________________________________________________________#
mus <- 0
sigmas <- c(seq(0.05, 0.45, by = 0.5),seq(0.5, 4.5, by = 0.5), seq(5,10,1))
n_samples <- 1e4
n_obs <- 50
# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(sigmas), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(sigmas), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(sigmas), ncol = 3),
                mhd_true = rep(0,length(sigmas)))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(sigmas), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- sapply(rep(sigmas[n],n_samples), function (x)
    rnorm(n_obs, mean = mus[1], sd = sigmas[1]), simplify = TRUE)
  # Calculat ethe MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[1]),sigmas[n]/sqrt(n_obs)))
  
  # Calculate the quantile of the MHD for each sample compared to pooled sample
  mhd_quantile_sample = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mean(x1)), sd(x1)/ sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(mu = mus[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
                         mhd_quantile_pop = mhd_quantile_pop, 
                         mhd_quantile_sample = mhd_quantile_sample)
  
  # Calculate summary stats
  df_sum$ci_95_sample[n,] = quantile(mhd_quantile_sample, probs = c(0.025,0.975))
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  df_sum$mhd_true[n] = sum(mhd > abs( mus[n]))
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[1])), n_samples, 
                                                    conf.level=0.95,methods="exact"))[c(4,5,6)]
  
}
df <- do.call("rbind", list_df)

# 0.95 quantile for each distribution
# MError of the 95% quantile 
p_b1 <- ggplot(df_sum, aes(x = sigma, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(mu)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_b1
save_plot(paste("figure/", fig_basename, "b t1e_MHD_over_sigma_mu_0.tiff", 
                sep = ""), p_b1, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)




# Visualize T1E (mu<MHD) under the alt hypothesis (mu=1), sigma swept
#
#______________________________________________________________________________#
mus <- 1
sigmas <- c(seq(0.05, 0.25, by = 0.05),seq(0.5, 4.5, by = 0.5), seq(5,10,1))
n_samples <- 1e4
n_obs <- 50
# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(sigmas), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(sigmas), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(sigmas), ncol = 3),
                mhd_true = rep(0,length(sigmas)))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(sigmas), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- sapply(rep(sigmas[n],n_samples), function (x)
    rnorm(n_obs, mean = mus[1], sd = sigmas[1]), simplify = TRUE)
  # Calculat ethe MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[1]),sigmas[n]/sqrt(n_obs)))
  
  # Calculate the quantile of the MHD for each sample compared to pooled sample
  mhd_quantile_sample = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mean(x1)), sd(x1)/ sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(mu = mus[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
                         mhd_quantile_pop = mhd_quantile_pop, 
                         mhd_quantile_sample = mhd_quantile_sample)
  
  # Calculate summary stats
  df_sum$ci_95_sample[n,] = quantile(mhd_quantile_sample, probs = c(0.025,0.975))
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  df_sum$mhd_true[n] = sum(mhd > abs( mus[1]))
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[1])), n_samples, 
                                                    conf.level=0.95,methods="exact"))[c(4,5,6)]
  
}
df <- do.call("rbind", list_df)

# 0.95 quantile for each distribution
# MError of the 95% quantile 
p_c1 <- ggplot(df_sum, aes(x = sigma, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(mu)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_c1
save_plot(paste("figure/", fig_basename, "c1 t1e_MHD_over_sigma_mu_1.tiff", 
                sep = ""), p_c1, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


        
        