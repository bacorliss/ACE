


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
                               mu_ov_sigmas = NULL, overrride = TRUE, rand.seed=0) {
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
  set.seed(rand.seed)
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



row_locate_binary_bounds <- function (xb){
  
  # Helper function to flip matrix left to right
  fliplr <- function(x) x[,ncol(x):1]
  
  # Code assumes that TRUE is on the left portion of matrix, if FALSE is there 
  # instead then invert
  if (sum(xb[,1]) < dim(xb)[1]) {
    xb = !xb
  }
  
  
  x <- matrix(as.numeric(xb),ncol = ncol(xb), nrow = nrow(xb))
  a<-t(apply(x,1,cumsum))
  b<-fliplr(t(apply(fliplr(!x),1,cumsum)))
  
  
  c <- b+a
  
  
  row_max_inds <- apply(c,1, function(x) mean(which(x == max(x))))
  return(row_max_inds)
  
}


locate_bidir_binary_thresh <- function(mus = NULL, sigmas, n_obs,temp_name, mu_ov_sigmas = NULL, rand.seed=0){
  #' Locate all 4 error transition boundaries
  #' Positive and negative direction, 0 and alpha error
  set.seed(rand.seed)
  
  n_cols <- max(c(length(mus), length(mu_ov_sigmas)))
  if (!is.null(mus)) {
    col_values = mus; varname <- "critical_mu"; neg_mu_ov_sigmas = mu_ov_sigmas;
  } else {
    col_values = mu_ov_sigmas; varname <- "critical_mu_over_sigma"; neg_mu_ov_sigmas = -mu_ov_sigmas;
  }
  
  # Run simulations calculating error of mmd with mu and sigma swept
  df_right <- stats_param_sweep(mus = mus, sigmas = sigmas, 
                                n_samples, n_obs, paste("temp/right_", temp_name,".rds"),
                                mu_ov_sigmas) 
  df_right$side <- as.factor("Right")
  
  # Run simulations calculating error of mmd with mu and sigma swept
  df_left <- stats_param_sweep(mus = mus, sigmas = sigmas, n_samples, n_obs, 
                               paste("temp/right_", temp_name,".rds"), 
                               mu_ov_sigmas= neg_mu_ov_sigmas) 
  df_left$side <- as.factor("Left")
  
  
  # Equivalence test versus middle column of same row
  p_threshold = 0.05 /(n_cols*length(sigmas))
  df_zero <- rbind(
    tibble(er = "0", side = "right", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_right$p_val_mmd_eq_zero  < p_threshold)),
    tibble(er = "alpha", side = "right", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_right$p_val_mmd_eq_alpha < p_threshold)))
  # Convert index of error rate transition to mu/sigma value
  df_zero[[varname]] <- 
    approx(x=1:n_cols,y = abs(col_values),
           xout = df_zero$critical_val_ind, 
           n = length(n_cols)*2L-1L)$y
  
  df_alpha <- rbind(
    tibble(er = "0", side = "left", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_left$p_val_mmd_eq_zero  < p_threshold)),
    tibble(er = "alpha", side = "left", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_left$p_val_mmd_eq_alpha  < p_threshold)))
  # Convert index of error rate transition to mu/sigma value
  df_alpha[[varname]] <- 
    approx(x=1:n_cols,y = abs(col_values),
           xout = df_alpha$critical_val_ind, 
           n = length(n_cols)*2L-1L)$y
  
  # Concatenate right and left dataframes of results
  df_crit <- rbind(df_zero,df_alpha)
  
  df_crit[[varname]][df_crit$side=="left"] = -
    abs(df_crit[[varname]][df_crit$side=="left"])
  # * abs() added just in case we are running the code multiple times when debugging
  
  
  return(df_crit)
  
}




rowcol_slopes <- function(m, row_vals, col_vals) {
  #' Calculate pearson correlation on a heatmap row by row and column by column,
  #' return bonferroni corrected p values
  #' 
  bonf_corr = dim(m)[1] + dim(m)[2]
  
  df_row <- data.frame(row_vals = row_vals, pearson_p = rep(0,dim(m)[1]), slope = rep(0,dim(m)[1]), 
                       slope_95L = rep(0,dim(m)[1]), slope_95U = rep(0,dim(m)[1]))
  
  for (r in seq( dim(m)[1])) {
    # Pearson p value
    df_row$pearson_p[r] = cor.test(unname(m[r,]),abs(col_vals),method = "pearson")$p.value/bonf_corr
    # Linear regression
    fit <- lm( y~x, data = data.frame(x=abs(col_vals), y = m[r,]))
    conf <- confint(fit, "x", level=1-0.05/bonf_corr)
    df_row$slope[r] <- fit$coefficients[2]
    df_row$slope_95L[r] <- conf[1]
    df_row$slope_95U[r] <- conf[2]
  }
  # df_row$sig_labels <- ifelse(row_pearson_p<0.05, "*","")
  
  
  df_col <- data.frame(col_vals = col_vals, pearson_p = rep(0,dim(m)[2]), slope = rep(0,dim(m)[2]), 
                       slope_95L = rep(0,dim(m)[2]), slope_95U = rep(0,dim(m)[2]))
  for (c in seq(dim(m)[2])){
    # Pearson p value
    df_col$pearson_p[c] = cor.test(unname(m[,c]),row_vals,method = "pearson")$p.value/bonf_corr
    # Linear regression
    fit <- lm( y~x, data = data.frame(x=abs(row_vals), y = m[,c]))
    conf <- confint(fit, "x", level=1-0.05/bonf_corr)
    df_col$slope[c] <- fit$coefficients[2]
    df_col$slope_95L[c]  <- conf[1]
    df_col$slope_95U[c]  <- conf[2]
  } 
  # df_col$col_sig_labels <- ifelse(col_pearson_p<0.05, "*","")
  
  slope_stats <- list(df_row, df_col)
  names(slope_stats)<-c("df_row", "df_col")
  return(slope_stats)
}


