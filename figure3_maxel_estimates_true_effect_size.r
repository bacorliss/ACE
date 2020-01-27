



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

source("R/mhd.r")


fig_basename = "f_3"
n_samples <- 1e3
rand_seed <-0


# 2D visualization of MHD differecne and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
# n_samples <- 1e4
n_obs <- 50

# Matrix diff and error of MHD
df_2d <- tibble(mean_diff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = length(mus)),
                mean_rdiff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = length(mus)),
                mean_mhd_error = matrix(0L, nrow = length(sigmas), ncol = length(mus)))

# Store results to disk since calculations are significant
set.seed(rand_seed)
if (!file.exists("temp/MHD_Error_2D.rds")) {
  
  # Generate random samples based on random mu values
  for (c in seq(1, length(mus), by = 1)) {
    for (r in seq(1, length(sigmas), by = 1)) {
      # Get samples, where each row is a seperate sample, columns are observations
      x1 <-
        matrix(rnorm(n_samples * n_obs, mus[c], sigmas[r]), nrow = n_obs)
      
      # Calculate the MHD
      mhd = apply(x1, 2, function (x)
        mhd_1sample_normal(x, alpha = 0.05))
      
      # Difference MHD to mu
      df_2d$mean_diff_mhd_mu[r, c] <- mean(mhd) - abs(mus[c])
      # Relative difference MHD to mu
      df_2d$mean_rdiff_mhd_mu[r, c] <- df_2d$mean_diff_mhd_mu[r, c] / sigmas[r]
      # Error rate MHD and mu
      df_2d$mean_mhd_error[r, c] <- sum(mhd < abs(mus[c])) / n_samples
      
    }
  }
  
  # Replace INF with NaNs
  df_2d$mean_rdiff_mhd_mu[!is.finite(df_2d$mean_rdiff_mhd_mu)] <- NaN
  # Save an object to a file
  saveRDS(df_2d, file = "temp/MHD_Error_2D.rds")
  
} else {
  # Restore the object
  df_2d <- readRDS(file = "temp/MHD_Error_2D.rds")
}


# https://sebastianraschka.com/Articles/heatmaps_in_r.html





# Visualize difference between MHD and mu
png(paste("figure/", fig_basename, "a1 MHD diff.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white","blue", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, .4, length=100), seq(0.41, 0.9, length=100),seq(0.91,2,length=100))             

heatmap.2(df_2d$mean_diff_mhd_mu, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow=rep("",length(sigmas)), labCol= round(mus,1),
          col = my_palette,
          density.info='none', scale="none",
          cexRow = 1, cexCol = 1,
          denscol="black", keysize=1,
          key.par=list(mar=c(3.5,0,3,0)),
          breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          key.title = expression(paste("Difference",~MHD~-abs(phantom(.)*mu*phantom(.)))),
          na.color = "black",
          key.xlab = "",
          main = NULL,
          xlab(expression(mu)),
          ylab(expression(sigma)))
dev.off()




# Visualize relative difference between MHD and mu
png(paste("figure/", fig_basename, "a2 MHD rdiff.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "red"))(n = 199)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, .3, length=100), seq(0.31, .4, length=100))             

heatmap.2(df_2d$mean_rdiff_mhd_mu, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow = rep("",length(sigmas)), labCol= round(mus,1),
          col = my_palette,
          density.info = 'none', scale="none",
          cexRow = 1, cexCol = 1,
          denscol = "black", keysize=1,
          key.par = list(mar=c(3.5,0,3,0)),
          breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          key.title = expression(paste("Relative Difference",~(MHD~-abs(~mu))/~sigma)),
          na.color = "black",
          key.xlab = "",
          main = NULL,
          xlab(expression(mu)),
          ylab(expression(sigma)))
dev.off()



# Visualize error rate of MHD
png(paste("figure/", fig_basename, "a3 MHD error rate 2D.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, 0.039, length=100), seq(0.04, 0.06, length=100),seq(0.061, 0.08,length=100)) 

color_cull <- function(x) x[seq(1,round(length(x)*.73), by = 1)]
heatmap.2(df_2d$mean_mhd_error, Rowv = FALSE, Colv = FALSE, trace = "none", dendrogram = "none", 
          labRow = rev(round(sigmas,1)), labCol = round(mus,1),
          col = my_palette,
          density.info='none', scale="none",
          cexRow = 1, cexCol = 1, denscol="black", keysize=1, 
          key.par = list( mar = c(3.5,0,3,0)),
          lmat = rbind(c(5, 4, 2), c(6, 1, 3)), margins = c(3,0),
          lhei = c(2, 6), lwid = c(1, 10, 1),
          breaks = col_breaks,
          key.title = "",
          key.xlab = "",
          main = NULL, 
          xlab( expression(mu)), 
          ylab( expression(sigma)))
dev.off()





# Visualize T1E (mu > MHD) under the null and alt hypothesis, mu swept (sigma=1)
#                                                                              #
#______________________________________________________________________________#

mus <- seq(-10, 10, by = .1)
sigmas <- c(1)
# n_samples <- 1e4
n_obs <- 50
# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(mus), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(mus), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(mus), ncol = 3),
                mhd_true = rep(0,length(mus)),
                diff_mhd_mu = matrix(0L, nrow = length(mus), ncol = 3))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(mus), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- matrix( rnorm( n_samples*n_obs, mus[n], sigmas[1]), nrow=n_obs )
  
  # Calculate the MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[n]),sigmas[1]/sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(mu = mus[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
           mhd_quantile_pop = mhd_quantile_pop, 
           mhd_quantile_sample = mhd_quantile_sample)
  
  # Calculate summary stats
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  df_sum$mhd_true[n] = sum(mhd>abs(mus[n]))
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[n])), n_samples, 
                                       conf.level=0.95,methods="exact"))[c(4,5,6)]
  # Record mean diff of MHD and mu and the 95% quantiles
  df_sum$diff_mhd_mu[n,] = c(mean(mhd - abs(mus[n])),quantile(mhd-abs(mus[n]),c(0.025,0.975)) )
  
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

# 95% quantile of the error for each distribution
p_a2 <- ggplot(df_sum, aes(x = mu, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(mu)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_a2
save_plot(paste("figure/", fig_basename, "a2 t1e_MHD_over_mu.tiff", 
                sep = ""), p_a2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# zoomed in 95% quantile of the error for each distribution
p_a3 <- ggplot(subset(df_sum,mu>=-0.5 & mu<=0.5), aes(x = mu, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha=.2) +
  xlab(expression(mu)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_a3
save_plot(paste("figure/", fig_basename, "a3 t1e_MHD_over_mu_zoomed.tiff", 
                sep = ""), p_a3, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# Difference between MHD and Mu 
p_a4 <- ggplot(df_sum, aes(x = mu, y = diff_mhd_mu[,1])) + 
  geom_hline(aes(yintercept = 0.0), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = diff_mhd_mu[,2], ymax = diff_mhd_mu[,3]),alpha = .2) +
  xlab(expression(mu)) + ylab(expression(M~H~D-~mu)) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_a4
save_plot(paste("figure/", fig_basename, "a4 diff_MHD_mu_over_mu_sigma_1.tiff", 
                sep = ""), p_a4, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)




# Visualize T1E (mu<MHD) under the null hypothesis (mu=0), sigma swept
#                                                                              #
#______________________________________________________________________________#

mus <- 0
sigmas <- seq(0.1, 2.5, by = 0.05)
# n_samples <- 1e4
n_obs <- 50
# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(sigmas), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(sigmas), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(sigmas), ncol = 3),
                mhd_true = rep(0,length(sigmas)),
                diff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = 3))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(sigmas), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- matrix( rnorm( n_samples*n_obs, mus[1], sigmas[n]), nrow=n_obs )
  
  # Calculat ethe MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[1]),sigmas[n]/sqrt(n_obs)))
  
  # Calculate the quantile of the MHD for each sample compared to pooled sample
  # mhd_quantile_sample = 
  #   sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mean(x1)), sd(x1)/ sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(sigma = sigmas[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
                         mhd_quantile_pop = mhd_quantile_pop)
  
  # Calculate summary stats
  #df_sum$ci_95_sample[n,] = quantile(mhd_quantile_sample, probs = c(0.025,0.975))
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  df_sum$mhd_true[n] = sum(mhd > abs( mus[1]))
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[1])), n_samples, 
                                                    conf.level=0.95,methods="exact"))[c(4,5,6)]
  # Record mean diff of MHD and mu and the 95% quantiles
  df_sum$diff_mhd_mu[n,] = c(mean(mhd - abs(mus[1])),quantile(mhd-abs(mus[1]),c(0.025,0.975)) )
}
df <- do.call("rbind", list_df)

# 0.95 quantile for each distribution
p_b1 <- ggplot(df_sum, aes(x = sigma, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(sigma)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_b1
save_plot(paste("figure/", fig_basename, "b t1e_MHD_over_sigma_mu_0.tiff", 
                sep = ""), p_b1, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

# Difference between MHD and Mu 
p_b2 <- ggplot(df_sum, aes(x = sigma, y = diff_mhd_mu[,1])) + 
  geom_hline(aes(yintercept = 0.0), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = diff_mhd_mu[,2], ymax = diff_mhd_mu[,3]),alpha = .2) +
  xlab(expression(sigma)) + ylab(expression(M~H~D-~mu)) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_b2
save_plot(paste("figure/", fig_basename, "b2 diff_MHD_mu_over_sigma_mu_0.tiff", 
                sep = ""), p_b2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# Difference between MHD and Mu 
p_b3 <- ggplot(df_sum, aes(x = sigma, y = diff_mhd_mu[,1]/sigma)) + 
  geom_hline(aes(yintercept = 0.0), color="red")+ geom_line() + 
  geom_ribbon(aes(ymin = diff_mhd_mu[,2]/sigma, ymax = diff_mhd_mu[,3]/sigma),alpha = .2) +
  xlab(expression(sigma)) + ylab(expression((M~H~D-~mu)/~sigma)) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_b3
save_plot(paste("figure/", fig_basename, "b3 rel_diff_MHD_mu_over_sigma_mu_0.tiff", 
                sep = ""), p_b3, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)



# Visualize T1E (mu<MHD) under the alt hypothesis (mu=1), sigma swept
#
#______________________________________________________________________________#
mus <- 1
sigmas <- c(seq(0.05, 0.25, by = 0.05),seq(0.5, 4.5, by = 0.5), seq(5,10,1))
n_obs <- 50
# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(sigmas), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(sigmas), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(sigmas), ncol = 3),
                mhd_true = rep(0,length(sigmas)),
                diff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = 3))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(sigmas), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- matrix( rnorm( n_samples*n_obs, mus[1], sigmas[n]), nrow=n_obs )

  # Calculat ethe MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[1]),sigmas[n]/sqrt(n_obs)))
  
  # Calculate the quantile of the MHD for each sample compared to pooled sample
  # mhd_quantile_sample = 
  #   sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mean(x1)), sd(x1)/ sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(mu = mus[1], sigma = sigmas[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
                         mhd_quantile_pop = mhd_quantile_pop, 
                         mhd_quantile_sample = mhd_quantile_sample)
  
  # Calculate summary stats
  # 95% CI of the MHD estimated from observations pooled across all samples
  #df_sum$ci_95_sample[n,] = quantile(mhd_quantile_sample, probs = c(0.025,0.975))
  # 95% CI of the MHD estimated from population parameters
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  # Binary is MHD "True" (Does it contain mu?)
  df_sum$mhd_true[n] = sum(mhd > abs( mus[1]))
  # Estimate confidence interval of binomial porportion of true MHDs
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[1])), n_samples, 
                                                    conf.level=0.95,methods="exact"))[c(4,5,6)]
  # Record mean diff of MHD and mu and the 95% quantiles
  df_sum$diff_mhd_mu[n,] = c(mean(mhd-mus[1]),quantile(mhd-mus[1],c(0.025,0.975)) )
}
df <- do.call("rbind", list_df)


# 95% quantile of T1 error rate
p_c1 <- ggplot(df_sum, aes(x = sigma, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(sigma)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_c1
save_plot(paste("figure/", fig_basename, "c1 t1e_MHD_over_sigma_mu_1.tiff", 
                sep = ""), p_c1, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# Difference between MHD and Mu 
p_c2 <- ggplot(df_sum, aes(x = sigma, y = diff_mhd_mu[,1])) + 
  geom_hline(aes(yintercept = 0.0), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = diff_mhd_mu[,2], ymax = diff_mhd_mu[,3]),alpha = .2) +
  xlab(expression(sigma)) + ylab(expression(M~H~D-~mu)) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_c2
save_plot(paste("figure/", fig_basename, "c2 diff_MHD_mu_over_sigma_mu_1.tiff", 
                sep = ""), p_c2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

# Difference between MHD and Mu 
p_c2 <- ggplot(df_sum, aes(x = sigma, y = diff_mhd_mu[,1]/sigma)) + 
  geom_hline(aes(yintercept = 0.0), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = diff_mhd_mu[,2]/sigma, ymax = diff_mhd_mu[,3]/sigma),alpha = .2) +
  xlab(expression(sigma)) + ylab(expression((M~H~D-~mu)/~sigma)) +
  theme_classic(base_size = 8) + theme(legend.position="none")
p_c2
save_plot(paste("figure/", fig_basename, "c2 rel_diff_MHD_mu_over_sigma_mu_1.tiff", 
                sep = ""), p_c2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)




        

# Test for distinct cut-off in T1E for mu/sigma
#
#______________________________________________________________________________#
n_samples <- 1e4
# Randomly generate large number of mus and sigmas (n_samples # of them)
mus <- 1
sigmas <- seq(1,10,0.1)
n_obs <- 50

# Matrix summarizing output
df_sum = tibble(mu = mus, sigma = sigmas, 
                ci_95_sample= matrix(0L, nrow = length(sigmas), ncol = 2),
                ci_95_pop= matrix(0L, nrow = length(sigmas), ncol = 2),
                mhd_true_ci = matrix(0L, nrow = length(sigmas), ncol = 3),
                mhd_true = rep(0,length(sigmas)),
                diff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = 3))

# Generate random samples based on random mu values
set.seed(rand_seed)
list_df <- list()
for (n in seq(1, length(sigmas), by=1)) {
  # Get samples, where each row is a seperate sample, columns are observations
  x1 <- matrix( rnorm( n_samples*n_obs, mus[1], sigmas[n]), nrow=n_obs )

  # Calculate the MHD
  mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha=0.05))
  
  # Calculate the quantile of the MHD for each sample compared to population
  mhd_quantile_pop = 
    sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mus[1]),sigmas[n]/sqrt(n_obs)))
  
  # # Calculate the quantile of the MHD for each sample compared to pooled sample
  # mhd_quantile_sample = 
  #   sapply(1:n_samples, function(i) pnorm(mhd[i], abs(mean(x1)), sd(x1)/ sqrt(n_obs)))
  
  # Compile output data
  list_df[[n]] <- tibble(mu = mus[1], sigma = sigmas[n],  x_bar = rowMeans(t(x1)), mhd = mhd,
                         mhd_quantile_pop = mhd_quantile_pop, 
                         mhd_quantile_sample = mhd_quantile_sample)
  
  # Calculate summary stats
  # 95% CI of the MHD estimated from observations pooled across all samples
  #df_sum$ci_95_sample[n,] = quantile(mhd_quantile_sample, probs = c(0.025,0.975))
  # 95% CI of the MHD estimated from population parameters
  df_sum$ci_95_pop[n,] = quantile(mhd_quantile_pop, probs = c(0.025,0.975))
  # Binary is MHD "True" (Does it contain mu?)
  df_sum$mhd_true[n] = sum(mhd > abs( mus[1]))
  # Estimate confidence interval of binomial porportion of true MHDs
  df_sum$mhd_true_ci[n,] = as.numeric(binom.confint(sum(mhd>abs(mus[1])), n_samples, 
                                                    conf.level=0.95,methods="exact"))[c(4,5,6)]
  # Record mean diff of MHD and mu and the 95% quantiles
  df_sum$diff_mhd_mu[n,] = c(mean(mhd-mus[1]),quantile(mhd-mus[1],c(0.025,0.975)) )
}
df <- do.call("rbind", list_df)




# Plot 95% quantile of MHD-mu versus mu/sigma
p_d1 <- ggplot(df_sum, aes(x = mu/sigma, y = diff_mhd_mu[,1])) + 
  geom_hline(aes(yintercept = 0.0), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = diff_mhd_mu[,2], ymax = diff_mhd_mu[,3]),alpha = .2) +
  xlab(expression(mu/sigma)) + ylab(expression(M~H~D-~mu)) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_d1
save_plot(paste("figure/", fig_basename, "d1 diff_MHD_mu_over_sigma_mu_1.tiff", 
                sep = ""), p_d1, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)



# Plot 95% quantile of T1 Error Rate versus mu/sigma
p_d2 <- ggplot(df_sum, aes(x = mu/sigma, y = 1-mhd_true_ci[,1])) + 
  geom_hline(aes(yintercept = 0.05), color="red") + geom_line() + 
  geom_ribbon(aes(ymin = 1-mhd_true_ci[,2], ymax = 1-mhd_true_ci[,3]),alpha = .2) +
  xlab(expression(sigma)) + ylab("MHD Type 1 Error") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_d2
save_plot(paste("figure/", fig_basename, "d2 t1e_MHD_over_sigma_mu_1.tiff", 
                sep = ""), p_d2, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)








# plot histograms to show range of samples for mu and sigma, and mu/sigma


