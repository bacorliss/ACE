

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
library(tidyr)

# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
library(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
library(TOSTER)
library(docstring)

# Helper Functions
source("R/mmd.r")
RowVar <- function(x, ...) rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)

# Test if a value is within an interval (inclusive)
within_ci <- function(ci, x) x >= ci[1] & x <= ci[2]

error_test_codes <-function(is_error_rate_zero, is_error_rate_alpha) {
#' Assign codes for error rate whether null hypothesis is rejected
#' 
#' @description Test if proportion of mmds that are above mu (error rate) are equal to 0.05 
#' and 0.00, return code do delineate combinations of both results
#'  (0) E = neither,  (1) E == 0.00
#'  (2) E == 0.05,  (3) E == both
  assign_code <- function(z,a) { if (z & a) {value=3} else if(z & !a) {value = 2}
    else if(!z & a) {value = 1} else {value = 0}; return(value)
  }
  error_rate_codes <- matrix(mapply(assign_code, is_error_rate_zero, is_error_rate_alpha, SIMPLIFY = TRUE,
         USE.NAMES = FALSE), nrow = dim(is_error_rate_zero)[1], ncol = dim(is_error_rate_zero)[2])
  colnames(error_rate_codes) <- colnames(is_error_rate_zero)
  rownames(error_rate_codes) <- rownames(is_error_rate_zero)
  return(error_rate_codes)
}

stats_param_sweep <- function( mus, sigmas, n_samples, n_obs, out_path, 
                              mu_ov_sigmas = NULL, overrride = TRUE) {
  #' Perform parameter sweep with specified mus and sigmas
  #' 
  #' @description QUantifies stats of a series of simulations of normal random 
  #' samples (n_samples) with a specified number of observations (n_obs) with 
  #' the specified values for mu (mus) and sigma (sigmas) for a normal 
  #' distribution.
  #'  
  #'  @param mus vector of mu values
  #'  @param sigmas vector of sigma values
  #'  @param n_samples number of samples (collection of measurements to simulate
  #'   one experiment)
  #'  @param n_obs number of observations
  #'  @param out_path output path to store results of calculation so simulations
  #'   do not have to be rerun
  #'  @param mu_ov_sigmas 
  #'  
  #'  @return dataframe that stores input parameters 

  # Store results to disk since calculations are significant
  set.seed(rand_seed)
  if (!file.exists(out_path) || overrride) {
    n_mus = max(c(length(mus), length(mu_ov_sigmas)))
    # Matrix diff and error of mmd
    
    if (!is.null(mus)) {dimnames = list(sigmas,mus)
    } else {dimnames = list(sigmas,mu_ov_sigmas)}

    df <- tibble(
      sigma = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mu = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mu_over_sigma = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                             dimnames = dimnames),
      mean_diff_mmd_mu = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                dimnames = dimnames),
      mean_rdiff_mmd_mu = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                 dimnames = dimnames),
      mean_mmd_error = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                              dimnames = dimnames),
      error_rate_tests = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                dimnames = dimnames),
      p_val_mmd_eq_zero = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                 dimnames = dimnames),
      p_val_mmd_eq_alpha = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                  dimnames = dimnames))
    # browser()
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
        # Calculate the mmd
        mmd = apply(x1, 2, function (x) mmd_normal_zdist(x, conf.level = 0.95) )
        
        # Difference mmd to mu
        df$mean_diff_mmd_mu[r, c] <- mean(mmd) - abs(mus[c])
        # Relative difference mmd to mu
        df$mean_rdiff_mmd_mu[r, c] <- df$mean_diff_mmd_mu[r, c] / sigmas[r]
        # Error rate mmd and mu
        n_successes <- sum(mmd < abs(mus[c]))
        df$mean_mmd_error[r, c] <- n_successes / n_samples
        # Caluate p-values for test against error rate == 0
        df$p_val_mmd_eq_zero[r, c] <- binom.test(
          n_successes, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
        # Caluate p-values for test against error rate == alpha
        df$p_val_mmd_eq_alpha[r, c] <- binom.test(
          n_successes, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
        # # Fraction of monte carlo tries where mmd was in error
        # df$mean_mmd_error[r, c] <- n_successes / n_samples
      }
    }
    # Replace INF with NaNs
    df$mean_rdiff_mmd_mu[!is.finite(df$mean_rdiff_mmd_mu)] <- NaN
    # Save an object to a file
    saveRDS(df, file = out_path)
  } else {
    # Restore the object
    df <- readRDS(file = out_path)
  }
  return(df)
}


# General variables
fig_num = "3"
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)

n_samples <- 1e3
rand_seed <- 0



# 2D visualization of mmd difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
n_obs <- 50

# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma.rds") 


# Difference MMD from Mu
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_diff_mmd_mu)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
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
                       limits=c(0,2)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a1 mmd diff.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 


# Difference MMD from Mu
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rdiff_mmd_mu)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c("blue","white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a2 mmd rdiff.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 


# Difference MMD from Mu
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_mmd_error)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
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
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a3 mmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 




# Difference MMD from Mu
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$p_val_mmd_eq_zero > 0.05, 
  df_results$p_val_mmd_eq_alpha > 0.05))) %>% gather(mu, z, -sigma)
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
  scale_fill_manual(values=c("white", "blue","red","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a4 mmd error test.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600) 





# 2D visualization of mmd difference and error rate over sigma and mu/sigma
#                                                                              #
#______________________________________________________________________________#
sigmas <- seq(.1, 5, by = .1)
mu_ov_sigmas <- seq (-.5, .5, by=0.01)
n_obs <- 50
set.seed(rand_seed)
# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(NULL, sigmas, n_samples, n_obs, 
                           "temp/mmd_Error_2D_mu_over_sigma_vs_sigma.rds", mu_ov_sigmas) 
# Difference MMD from Mu
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$p_val_mmd_eq_zero > 0.05, 
  df_results$p_val_mmd_eq_alpha > 0.05))) %>% gather(mu, z, -sigma)
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
  scale_fill_manual(values=c("white", "blue","red","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a5 mmd error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600) 






# Identify location of coverage error boundaries with mu space
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
sigmas <- seq(.1, 5, by = .1); 
mus <- seq (.1/5, 2, by=0.01)

crit_df <- locate_bidir_binary_thresh(sigmas = sigmas, n_obs = n_obs, 
                                      temp_name = "mmd_Error_mu_over_sigma_vs_sigma.rds", 
                                      mu_ov_sigmas = mu_ov_sigmas)


# Box and Whiskers of Coverage Error Transition Region
p <- ggplot(crit_df, aes(x=er,  color = er, group = interaction(er, side), 
                         y = critical_mu_over_sigma)) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") + 
  ylab(expression(mu/sigma))
p
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_f mmd transition values.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(critical_mu_over_sigma ~ er + side, data = crit_df)
summary(res.aov2)



# Identify location of coverage error boundaries wiht mu/sigma
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
sigmas <- seq(.1, 5, by = .1); 
mu_ov_sigmas <- seq (0.15, 0.35, by=0.001)

crit_df <- locate_bidir_binary_thresh(sigmas = sigmas, n_obs = n_obs, 
                                      temp_name = "mmd_Error_mu_over_sigma_vs_sigma.rds", 
                                      mu_ov_sigmas = mu_ov_sigmas)


# Box and Whiskers of Coverage Error Transition Region
p <- ggplot(crit_df, aes(x=er,  color = er, group = interaction(er, side), 
                         y = critical_mu_over_sigma)) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") + 
  ylab(expression(mu/sigma))
p
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_f mmd transition values.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(critical_mu_over_sigma ~ er + side, data = crit_df)
summary(res.aov2)







locate_bidir_binary_thresh <- function(sigmas,n_obs,temp_name, mu_ov_sigmas){
  #' Locate all 4 error transition boundaries
  #' Positive and negative direction, 0 and alpha error
  
  # Run simulations calculating error of mmd with mu and sigma swept
  right_mu_ov_sigmas <- mu_ov_sigmas
  df_right <- stats_param_sweep(
    NULL, sigmas, n_samples, n_obs, paste("temp/right_", temp_name,".rds"),
    right_mu_ov_sigmas) 
  df_right$side <- as.factor("Right")

  # Run simulations calculating error of mmd with mu and sigma swept
  left_mu_ov_sigmas <- -mu_ov_sigmas
  df_left <- stats_param_sweep(
    NULL, sigmas, n_samples, n_obs, paste("temp/right_", temp_name,".rds"),
    left_mu_ov_sigmas) 
  df_left$side <- as.factor("Left")
  

  # Equivalence test versus middle column of same row
  p_threshold = 0.05 /length(right_mu_ov_sigmas)
  zero_df <- rbind(
    tibble(er = "0", side = "right", sigma = sigmas, 
           critical_mu_over_sigma_ind = apply(
             df_right$p_val_mmd_eq_zero  < p_threshold, 1, 
             locate_binary_thresh, vals = c(FALSE,TRUE))),
    tibble(er = "alpha", side = "right", sigma = sigmas, 
           critical_mu_over_sigma_ind = apply(
             df_right$p_val_mmd_eq_alpha  < p_threshold,1,
             locate_binary_thresh, vals = c(FALSE,TRUE))))
  # COnvert index of error rate transition to mu/sigma value
  zero_df$critical_mu_over_sigma <- 
    approx(x=seq_along(right_mu_ov_sigmas),y = abs(right_mu_ov_sigmas),
           xout = zero_df$critical_mu_over_sigma_ind, 
           n = length(mu_ov_sigmas)*2L-1L)$y
  
  alpha_df <- rbind(
    tibble(er = "0", side = "left", sigma = sigmas, 
           critical_mu_over_sigma_ind = apply(
             df_left$p_val_mmd_eq_zero  < p_threshold, 1, 
             locate_binary_thresh, vals = c(FALSE,TRUE))),
    tibble(er = "alpha", side = "left", sigma = sigmas, 
           critical_mu_over_sigma_ind = apply(
             df_left$p_val_mmd_eq_alpha  < p_threshold, 1, 
             locate_binary_thresh, vals = c(FALSE,TRUE))))
  # COnvert index of error rate transition to mu/sigma value
  alpha_df$critical_mu_over_sigma <- 
    approx(x=seq_along(left_mu_ov_sigmas),y = abs(left_mu_ov_sigmas),
           xout = alpha_df$critical_mu_over_sigma_ind, 
           n = length(mu_ov_sigmas)*2L-1L)$y
  # Concatenate right and left dataframes of results
  crit_df <- rbind(zero_df,alpha_df)
  
  return(crit_df)
  
}


locate_binary_thresh <- function(vect, vals) {
  #' For a binary vector, find a threshold that separates TRUE and FALSE with balanced
  #' degree of false positives.
  #' @param vect is input vectors of TRUE and FALSE
  #' vals: TRUe or FALSE for each side of vector
  #' Example input
  #' vect = df_right$p_val_mmd_eq_zero[1,] < 0.05
  #' vals = c(FALSE,TRUE)
  
  vector_index = seq_along(vect)
  # score each position based on  1/x distance from threshold
  right_score = rep(NaN,length(vector_index))
  left_score = rep(NaN,length(vector_index))
  score = rep(NaN,length(vector_index))
  
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

# Take test 0 and test 0.05 matrices and splits them in middle, locate transition boundaries
# Location: +mu, - mu, Error:0, 95









