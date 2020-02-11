

# p.adjust(p, method = p.adjust.methods, n = length(p))


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
# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
library(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
library(TOSTER)


# Helper Functions
source("R/mhd.r")
RowVar <- function(x, ...) rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)

# Test if a value is within an interval (inclusive)
within_ci <- function(ci, x) x >= ci[1] & x <= ci[2]

error_test_codes <-function(is_error_rate_zero, is_error_rate_alpha) {
# Test if proportion of MHDs that are above mu (error rate) are equal to 0.05 
# and 0.00, return code do delineate combinations of both results
#  (0) E = neither,  (1) E == 0.00
#  (2) E == 0.05,  (3) E == both
  assign_code <- function(z,a) { if (z & a) {value=3} else if(z & !a) {value = 2}
    else if(!z & a) {value = 1} else {value = 0}; return(value)
  }
  error_rate_codes <- matrix(mapply(assign_code, is_error_rate_zero, is_error_rate_alpha, SIMPLIFY = TRUE,
         USE.NAMES = FALSE), nrow = dim(is_error_rate_zero)[1], ncol = dim(is_error_rate_zero)[2])
  return(error_rate_codes)
}

stats_param_sweep <- function(mus, sigmas, n_samples, n_obs, out_path, mu_ov_sigmas=NULL) {
  
  # Store results to disk since calculations are significant
  set.seed(rand_seed)
  if (!file.exists(out_path)) {
    n_mus = max(c(length(mus), length(mu_ov_sigmas)))
    # Matrix diff and error of MHD
    df <- tibble(
      sigma = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mu = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mu_over_sigma = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mean_diff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mean_rdiff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mean_mhd_error = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      error_rate_tests = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      p_val_mhd_eq_zero = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      p_val_mhd_eq_alpha = matrix(0L, nrow = length(sigmas), ncol = n_mus))
    
    # Generate random samples based on random mu values
    for (r in seq(1, length(sigmas), by = 1)) {
      # With a vector of mu/sigma and sigma known, calculate mu
      if (!is.null(mu_ov_sigmas)) { mus = mu_ov_sigmas * sigmas[r]}
      
      for (c in seq(1, length(mus), by = 1)) {
        
        df$sigma[r, c] = sigmas[r]
        df$mu[r, c] = mus[c]
        df$mu_over_sigma[r, c] = mus[c]/sigmas[r]
        
        # Get samples, where each row is a seperate sample, columns are observations
        x1 <- matrix(rnorm(n_samples * n_obs, mus[c], sigmas[r]), nrow = n_obs)
        # Calculate the MHD
        mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha = 0.05))
        
        # Difference MHD to mu
        df$mean_diff_mhd_mu[r, c] <- mean(mhd) - abs(mus[c])
        # Relative difference MHD to mu
        df$mean_rdiff_mhd_mu[r, c] <- df$mean_diff_mhd_mu[r, c] / sigmas[r]
        # Error rate MHD and mu
        n_successes <- sum(mhd < abs(mus[c]))
        df$mean_mhd_error[r, c] <- n_successes / n_samples
        # Caluate p-values for test against error rate == 0
        df$p_val_mhd_eq_zero[r, c] <- binom.test(
          n_successes, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
        # Caluate p-values for test against error rate == alpha
        df$p_val_mhd_eq_alpha[r, c] <- binom.test(
          n_successes, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
        # # Fraction of monte carlo tries where MHD was in error
        # df$mean_mhd_error[r, c] <- n_successes / n_samples
      }
    }
    # Replace INF with NaNs
    df$mean_rdiff_mhd_mu[!is.finite(df$mean_rdiff_mhd_mu)] <- NaN
    # Save an object to a file
    saveRDS(df, file = out_path)
  } else {
    # Restore the object
    df <- readRDS(file = out_path)
  }
  return(df)
}


# General variables
fig_basename = "f_3"
n_samples <- 1e3
rand_seed <- 0




# 2D visualization of MHD difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
# n_samples <- 1e4
n_obs <- 50

# Run simulations calculating error of MHD with mu and sigma swept
df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/MHD_Error_2D_mu_vs_sigma.rds") 

# Visualize difference between MHD and mu
# https://sebastianraschka.com/Articles/heatmaps_in_r.html
png(paste("figure/", fig_basename, "a1 MHD diff.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white","blue", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, .4, length=100), seq(0.41, 0.9, length=100),seq(0.91,2,length=100))             
heatmap.2(df_results$mean_diff_mhd_mu, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow=rep("",length(sigmas)), labCol= round(mus,1),
          col = my_palette, density.info='none', scale="none",
          cexRow = 1, cexCol = 1, denscol="black", keysize=1,
          key.par=list(mar=c(3.5,0,3,0)), breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          key.title = "", na.color = "black", key.xlab = "", main = NULL,
          xlab(expression(mu)), ylab(expression(sigma)))
dev.off()


# Visualize relative difference between MHD and mu
png(paste("figure/", fig_basename, "a2 MHD rdiff.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "red"))(n = 199)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, .3, length=100), seq(0.31, .4, length=100))             
heatmap.2(df_results$mean_rdiff_mhd_mu, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow = rep("",length(sigmas)), labCol= round(mus,1),
          col = my_palette, density.info = 'none', scale="none",
          cexRow = 1, cexCol = 1, denscol = "black", keysize=1,
          key.par = list(mar=c(3.5,0,3,0)), breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          key.title = "", na.color = "black", key.xlab = "", main = NULL,
          xlab(expression(mu)), ylab(expression(sigma)))
dev.off()



# Visualize error rate of MHD
png(paste("figure/", fig_basename, "a3 MHD error rate 2D.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, 0.039, length=100), seq(0.04, 0.06, length=100),seq(0.061, 0.08,length=100)) 
color_cull <- function(x) x[seq(1,round(length(x)*.73), by = 1)]
heatmap.2(df_results$mean_mhd_error, Rowv = FALSE, Colv = FALSE, trace = "none", dendrogram = "none", 
          labRow = rev(round(sigmas,1)), labCol = round(mus,1),
          col = my_palette, density.info='none', scale="none",
          cexRow = 1, cexCol = 1, denscol="black", keysize=1, 
          key.par = list( mar = c(3.5,0,3,0)),
          lmat = rbind(c(5, 4, 2), c(6, 1, 3)), margins = c(3,0),
          lhei = c(2, 6), lwid = c(1, 10, 1),
          breaks = col_breaks, key.title = "",key.xlab = "", main = NULL, 
          xlab( expression(mu)), ylab( expression(sigma)))
dev.off()



# 2D Heatmap of error (1) p == 0.05), (2) p == 1, (3) p == 0.05 & p == 1
png(paste("figure/", fig_basename, "a4 MHD error test.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white","blue", "red", "purple"))(n = 399)
col_breaks = c(seq(0, .5, length=100), seq(0.6, 1.5, length=100), 
               seq(1.6,2.5,length=100), seq(2.6,3.5,length=100))    
heatmap.2(error_test_codes(df_results$p_val_mhd_eq_zero > 0.05, df_results$p_val_mhd_eq_alpha > 0.05),
          Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow=rev(round(sigmas,1)), labCol= round(mus,1),
          col = my_palette, breaks = col_breaks, density.info='none', scale="none",
          cexRow = 1, cexCol = 1, denscol="black", keysize=1,
          key = FALSE, lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          # key.title = expression(paste("Difference",~MHD~-abs(phantom(.)*mu*phantom(.)))),
          na.color = "black", main = NULL,
          xlab(expression(mu)), ylab(expression(sigma)))
dev.off()





# 2D visualization of MHD difference and error rate over sigma and mu/sigma
#                                                                              #
#______________________________________________________________________________#
sigmas <- seq(.1, 5, by = .1)
mu_ov_sigmas <- seq (-.5, .5, by=0.01)
n_obs <- 50
set.seed(rand_seed)
# Run simulations calculating error of MHD with mu and sigma swept
df_results <- stats_param_sweep(NULL, sigmas, n_samples, n_obs, 
                           "temp/MHD_Error_2D_mu_over_sigma_vs_sigma.rds", mu_ov_sigmas) 

# 2D Heatmap of error (1) p == 0.05), (2) p == 1, (3) p == 0.05 & p == 1
png(paste("figure/", fig_basename, "a5 MHD error test mu_over_sigma.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white","blue", "red", "purple"))(n = 399)
col_breaks = c(seq(0, .5, length=100), seq(0.6, 1.5, length=100), 
               seq(1.6,2.5,length=100), seq(2.6,3.5,length=100))  
heatmap.2(
  error_test_codes(df_results$p_val_mhd_eq_zero > 0.05, df_results$p_val_mhd_eq_alpha > 0.05),
  Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
  labRow=rev(round(sigmas,1)), labCol= round(mu_ov_sigmas,2),
  col = my_palette, breaks = col_breaks,
  density.info = 'none', scale = "none", cexRow = 1, cexCol = 1,
  denscol = "black", keysize = 1, key = FALSE,
  lmat = rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
  lhei = c(2, 6), lwid = c(1, 10, 1),
  na.color = "black", main = NULL,
  xlab(expression(mu)), ylab(expression(sigma)))
dev.off()





# 2D visualization of MHD difference and error rate over sigma and mu/sigma
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
sigmas <- seq(.1, 5, by = .1); 

# Run simulations calculating error of MHD with mu and sigma swept
right_mu_ov_sigmas <- seq (0.15, 0.35, by=0.001)
df_right <- stats_param_sweep(
  NULL, sigmas, n_samples, n_obs, "temp/MHD_Error_right_mu_over_sigma_vs_sigma.rds", right_mu_ov_sigmas) 
df_right$side <- as.factor("Right")

# Run simulations calculating error of MHD with mu and sigma swept
left_mu_ov_sigmas <- seq (-0.15, -0.35, by=-0.001)
df_left <- stats_param_sweep(
  NULL, sigmas, n_samples, n_obs, "temp/MHD_Error_left_mu_over_sigma_vs_sigma.rds", left_mu_ov_sigmas) 
df_left$side <- as.factor("Left")

false_positive_min_threshold <- function(vect, vals) {
  # For a binary vector, find a threshold that seperates TRUE and FALSE with balanced
  # degree of false positives.
  # vect is input vectors of TRUE and FALSE
  # vals: TRUe or FALSE for each side of vector
  # Example input
  # vect = df_right$p_val_mhd_eq_zero[1,] < 0.05
  # vals = c(FALSE,TRUE)
  
  vector_index = seq_along(vect)
  # score each position based on  1/x distance from threshold
  
  for (i in seq(1, length(vect)-1)) {
    right_index = seq(i+1, length(vect),by = 1)
    # Get distances for all false positives to the right
    right_false_dist <- (right_index - i)[(vect == vals[2])[right_index]]
    # Compute score to the right of threshold
    right_score[i] <- sum(right_false_dist ^ 2)
    # Get distances for all false positives to the right
    left_index = seq(1, i, by = 1)
    left_false_dist <- (i - left_index)[(vect == vals[1])[left_index]]
    # Compute score to the right of threshold
    left_score[i] <- sum(left_false_dist ^ 2)
    score[i] <- abs(right_score[i] - left_score[i])
    
  }
  
  balance_ind = which.min(score)
  return(balance_ind)
}



 
# Equivalence test versus middle column of same row
p_threshold = 0.05 /length(right_mu_ov_sigmas)
zero_df <- rbind(
  tibble(er = "0", side = "right", sigma = sigmas, 
         critical_ind_mu_over_sigma = apply(
           df_right$p_val_mhd_eq_zero  < p_threshold, 1, 
           false_positive_min_threshold, vals = c(FALSE,TRUE))),
  tibble(er = "alpha", side = "right", sigma = sigmas, 
         critical_ind_mu_over_sigma = apply(
           df_right$p_val_mhd_eq_alpha  < p_threshold,1,
           false_positive_min_threshold, vals = c(FALSE,TRUE))))
zero_df$critical_mu_over_sigma <- 
  approx(x=seq_along(right_mu_ov_sigmas),y = abs(right_mu_ov_sigmas),
         xout = zero_df$critical_ind_mu_over_sigma, 
         n = length(mu_ov_sigmas)*2L-1L)$y

alpha_df <- rbind(
  tibble(er = "0", side = "left", sigma = sigmas, 
         critical_ind_mu_over_sigma = apply(
           df_left$p_val_mhd_eq_zero  < p_threshold, 1, 
           false_positive_min_threshold, vals = c(FALSE,TRUE))),
  tibble(er = "alpha", side = "left", sigma = sigmas, 
         critical_ind_mu_over_sigma = apply(
           df_left$p_val_mhd_eq_alpha  < p_threshold, 1, 
           false_positive_min_threshold, vals = c(FALSE,TRUE))))
alpha_df$critical_mu_over_sigma <- 
  approx(x=seq_along(left_mu_ov_sigmas),y = abs(left_mu_ov_sigmas),
         xout = alpha_df$critical_ind_mu_over_sigma, 
         n = length(mu_ov_sigmas)*2L-1L)$y
# Concatenate right and left datafgrames of results
crit_df <- rbind(zero_df,alpha_df)



# heatmap.2(
#   matrix(as.numeric(mc_eq_alpha), nrow = length(sigmas), ncol = length(mu_ov_sigmas)),
#   Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
#   labRow=rev(round(sigmas,1)), labCol= round(mu_ov_sigmas,2),
#   density.info = 'none', scale = "none", cexRow = 1, cexCol = 1,
#   denscol = "black", keysize = 1, key = FALSE,
#   lmat = rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
#   lhei = c(2, 6), lwid = c(1, 10, 1),
#   na.color = "black", main = NULL,
#   xlab(expression(mu)), ylab(expression(sigma)))
# dev.off()



# Basic violin plot
p <- ggplot(crit_df, aes(x=er,  color = er, group = interaction(er, side), 
                         y = critical_mu_over_sigma)) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") + 
  ylab(expression(mu/sigma))
p
save_plot(paste("figure/", fig_basename, "f MHD transition values.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(critical_mu_over_sigma ~ er + side, data = crit_df)
summary(res.aov2)