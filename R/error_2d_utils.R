
# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(foreach)
p_load(doParallel)
p_load(tidyr)

source("R/parallel_utils.R")
source("R/mmd.R")
source("R/ldm.R")
source("R/row_effect_sizes.R")
mmd_functions <- parse_functions_source("R/mmd.R")
ldm_functions <- parse_functions_source("R/ldm.R")
row_effect_size_functions <- parse_functions_source("R/row_effect_sizes.R")

RowVar <- function(x, ...) rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
# Test if a value is within an interval (inclusive)
within_ci <- function(ci, x) x >= ci[1] & x <= ci[2]




slopes_by_rowcol <- function(m, row_vals, col_vals) {
  #' Calculate pearson correlation on a 2d heatmap by row and by column,
  #' return bonferroni corrected p values to investigate trends. Row slope 
  #' calculation is done with  absolute value of the variable for rows because
  #' these heatmaps are symmetric about the origin.
  #' 
  #' @param m 2d matrix of heatmap
  #' @param row_vals values for each column of matrix
  #' @param col_cols values for each row of matrix
  #' @return slope_stats a 2 element names list that each contain a dataframe,
  #' one for the slopes for rows (df_row), and columns (df_col).
  bonf_corr = dim(m)[1] + dim(m)[2]
  
  df_row <- data.frame(row_vals = row_vals, pearson_p = rep(0,dim(m)[1]), slope = rep(0,dim(m)[1]), 
                       slope_95L = rep(0,dim(m)[1]), slope_95U = rep(0,dim(m)[1]))
  
  for (r in seq( dim(m)[1])) {
    # Pearson p value
    if (all(is.nan(unname(m[r,])))) {
      df_row$pearson_p[r] <- NaN
      df_row$slope_95L[r] <- NaN
      df_row$slope_95U[r] <- NaN
    } else {
      df_row$pearson_p[r] = cor.test(unname(m[r,]),abs(col_vals),method = "pearson")$p.value/bonf_corr
      # Linear regression
      fit <- lm( y~x, data = data.frame(x=abs(col_vals), y = m[r,]))
      conf <- confint(fit, "x", level=1-0.05/bonf_corr)
      df_row$slope[r] <- fit$coefficients[2]
      df_row$slope_95L[r] <- conf[1]
      df_row$slope_95U[r] <- conf[2]
    }
  }
  # df_row$sig_labels <- ifelse(row_pearson_p<0.05, "*","")
  
  
  df_col <- data.frame(col_vals = col_vals, pearson_p = rep(0,dim(m)[2]), slope = rep(0,dim(m)[2]), 
                       slope_95L = rep(0,dim(m)[2]), slope_95U = rep(0,dim(m)[2]))
  for (c in seq(dim(m)[2])){
    
    if (all(is.nan(unname(m[,c])))) {
      df_col$pearson_p[c] <- NaN
      df_col$slope_95L[c] <- NaN
      df_col$slope_95U[c] <- NaN
    } else {
      # Pearson p value
      df_col$pearson_p[c] = cor.test(unname(m[,c]),row_vals,method = "pearson")$p.value/bonf_corr
      # Linear regression
      fit <- lm( y~x, data = data.frame(x=abs(row_vals), y = m[,c]))
      conf <- confint(fit, "x", level=1-0.05/bonf_corr)
      df_col$slope[c] <- fit$coefficients[2]
      df_col$slope_95L[c]  <- conf[1]
      df_col$slope_95U[c]  <- conf[2]
      
    }
    
  } 
  # df_col$col_sig_labels <- ifelse(col_pearson_p<0.05, "*","")
  
  slope_stats <- list(df_row, df_col)
  names(slope_stats)<-c("df_row", "df_col")
  return(slope_stats)
}





error_test_codes <-function(is_error_rate_zero, is_error_rate_alpha) {
  #' @description Assign codes for error rate whether null hypothesis is rejected
  #' for use in a heatmap of hypothesis outcomes. Test if proportion of mmds that
  #' are above mu (error rate) are equal to 0.05 and 0.00, return code do delineate
  #' combinations of both results:
  #'  (0) E = neither,  (1) E == 0.00
  #'  (2) E == 0.05,  (3) E == both
  assign_code <- function(z,a) { if (z & a) {value=3} else if(z & !a) {value = 2}
    else if(!z & a) {value = 1} else {value = 0}; return(value)
  }
  error_rate_codes <- matrix(mapply(assign_code, is_error_rate_zero, 
                                    is_error_rate_alpha, SIMPLIFY = TRUE,
                                    USE.NAMES = FALSE), nrow = dim(is_error_rate_zero)[1], 
                             ncol = dim(is_error_rate_zero)[2])
  colnames(error_rate_codes) <- colnames(is_error_rate_zero)
  rownames(error_rate_codes) <- rownames(is_error_rate_zero)
  return(error_rate_codes)
}




quant_coverage_errors <-
  function( mus_ao, sigmas_ao, n_samples, n_obs, out_path, 
            mu_vsigmas_ao = NULL, overwrite = FALSE, rand.seed=0,
            mus_a = rep(0, max(c(length(mus_ao),length(mu_vsigmas_ao)))),
            sigmas_a = rep(0, max(c(length(mus_ao),length(mu_vsigmas_ao)))),
            is_parallel_proc = TRUE) {
    #' Perform parameter sweep with specified mus and sigmas
    #' 
    #' @description QUantifies stats of a series of simulations of normal random 
    #' samples (n_samples) with a specified number of observations (n_obs) with 
    #' the specified values for mu (mus_ao) and sigma (sigmas_ao) for a normal 
    #' distribution.
    #'  
    #'  @param mus_ao offset value between mus_a and mus_b, see project notes why 
    #'  this is used instead of mus_d. mus_b = mus_a + mus_ao
    #'  @param sigmas_ao offset between sigmas_a and sigmas_b. 
    #'  sigmas_b = sigmas_ao - sigmas_a
    #'  @param n_samples number of samples (collection of measurements to simulate
    #'   one experiment)
    #'  @param n_obs number of observations
    #'  @param out_path output path to store results of calculation so simulations
    #'   do not have to be rerun
    #'  @param mu_vsigmas_ao used in place of the mus_ao vector when looking at 
    #'  normalized values of mu
    #'  
    #'  @return df dataframe that stores population parameters and statistics
    
    dir.create(dirname(out_path),showWarnings = FALSE, recursive = TRUE)
    # browser();
    
    # save(list = ls(all.names = TRUE), file = "temp/quant_coverage_errors.RData",envir = environment())
    # load(file = "temp/quant_coverage_errors.RData")
    
    # Calculate number of simulations, which is either length of mus_ao or mu_vsigmas_ao
    n_mus_ao = max(c(length(mus_ao), length(mu_vsigmas_ao)))
    # Fill out vectors if single number specified
    if (length(mus_a)==1) {mus_a = rep(mus_a,n_mus_ao)}
    if (length(sigmas_a)==1) {sigmas_a = rep(sigmas_a,length(sigmas_ao))}

    # Store results to disk since calculations are significant
    set.seed(rand.seed)
    if (!file.exists(out_path) || overwrite) {

      # Matrix diff and error of mmd
      if (!is.null(mus_ao)) {dimnames = list(sigmas_ao,mus_ao)
      } else {dimnames = list(sigmas_ao,mu_vsigmas_ao)}
      
      # Initialize matrix in data form
      param_col_list <- c("mu_a", "sigma_a","mu_ao", "sigma_ao","mu_d", "sigma_d", 
                          "mu_vsigma_d","rmu","fract_neg_x_bar_a")
      init_mat <- matrix(NA, nrow = length(sigmas_ao), ncol = n_mus_ao, dimnames = dimnames)
      df_mat <- tibble(mu_a = init_mat);
      for (n in seq_along(param_col_list)) {df_mat[[param_col_list[n]]] <- init_mat }
      
      for (r in seq(1, length(sigmas_ao), by = 1)) {
        # Calculate mu if only mu/sigma and sigma are known
        if (!is.null(mu_vsigmas_ao)) { mus_ao = mu_vsigmas_ao * sigmas_ao[r]}
        for (c in seq(1, length(mus_ao), by = 1)) {
          df_mat$mu_ao[r,c]    <- mus_ao[c]
          df_mat$sigma_ao[r,c] <- sigmas_ao[r]
          # Control group
          df_mat$mu_a[r,c]     <- mus_a[c]
          df_mat$sigma_a[r,c]  <- sigmas_a[r]
          # Difference group
          df_mat$mu_d[r,c]     <- df_mat$mu_ao[r,c]
          df_mat$sigma_d[r,c]  <- sqrt(df_mat$sigma_a[r,c]^2 + (df_mat$sigma_a[r,c] +
                                                                  df_mat$sigma_ao[r,c])^2 )
          df_mat$mu_vsigma_d[r,c] <- df_mat$mu_d[r,c]/df_mat$sigma_d[r,c]
          df_mat$rmu[r,c] <- df_mat$mu_d[r,c] / df_mat$mu_a[r,c]
        }
      } 
      
      # Linearize matrix dataframe for parallel processing
      df_lin <- tibble(mu_a = as.vector(init_mat));
      for (n in seq_along(param_col_list)) { df_lin[[param_col_list[n]]] <- 
        as.vector(df_mat[[param_col_list[n]]]) }
      
      
      # Process parallel or serially
      if (is_parallel_proc) { 
        print("parallell")
        cl = makeCluster(detectCores()[1]-1)
        registerDoParallel(cl)
        df_lin2 <- foreach(n = seq(1,length(sigmas_ao)*length(mus_ao),1),
                           .export = c(mmd_functions, ldm_functions,row_effect_size_functions, 
                                       "quant_coverage_error","quant_error_rate"), 
                           .combine = rbind, .packages = c("tidyr")) %dopar% {
                             #calling a function
                             tempMatrix <- quant_coverage_error(df_lin[n,],n_samples, n_obs) 
                             tempMatrix
                           }
        stopCluster(cl)
      } else {            # Process effect sizes serially
        # browser();
        df_list <- list()
        for (n in seq(1,length(sigmas_ao)*length(mus_ao),1)) {
          print(n)
          df_list[[n]] <- quant_coverage_error(df_lin[n,],n_samples, n_obs) 
        }
        # save(list = ls(all.names = TRUE), file = "temp/quant_coverage_errors.RData",
        #      envir = environment())
        df_lin2 <- do.call("rbind", df_list)
      }
      
      # Convert linear data frame to dataframe of matrices
      df_mat2 <- df_mat; col_list <- colnames(df_lin2)
      for (n in seq_along(col_list)) { 
        df_mat2[[col_list[n]]] <- matrix(df_lin2[[col_list[n]]], nrow = length(sigmas_ao), 
                                         ncol = n_mus_ao, dimnames = dimnames)
      }
      
      
      # Save an object to a file
      saveRDS(df_mat2, file = out_path)
    } else {
      # Restore the object
      df_mat2 <- readRDS(file = out_path)
    }
    save(list = ls(all.names = TRUE), file = "temp/quant_coverage_errors.RData",
         envir = environment())
    return(df_mat2)
  }


quant_coverage_error <-  function(df, n_samples, n_obs) {
  #' @description Calculates the coverage error of x_bar, rx_bar, mmd, and rmmd
  #'
  #' @param df: a single row dataframe generated from agreement_contest
  #'  library, returned from generateExperiment_Data()
  #' @param n_samples: number of samples drawn to evaluate converage error
  #' @param n_obs: 
  #' @return null, exports figures to disk
  
  # save(list = ls(all.names = TRUE), file = "temp/quant_coverage_error.RData",envir = environment())
  # load(file = "temp/quant_coverage_error.RData")
 
  # browser();
  
  # Control group (For use with two sample cases)
  x_a <- matrix(rnorm(n_samples * n_obs, df$mu_a, df$sigma_a), ncol = n_obs)
  # Difference group (for simplicity experiment sample not calculated)
  x_d <- matrix(rnorm(n_samples * n_obs, df$mu_d, df$sigma_d), ncol = n_obs)
  # Each row is a separate sample, columns are observations
  
  ci_mean = row_ci_mean_2s_zdist(m1 = x_a, m2 = x_a+x_d)
  df_init = data.frame(xbar_dm = rowMeans(x_d), mmd = row_mmd_2s_zdist(x_a, x_a+x_d), 
                    ldm = row_ldm_2s_zdist(x_a, x_a+x_d), ci_lower_z = ci_mean$ci_lower, 
                    ci_upper_z = ci_mean$ci_upper, mu_dm = df$mu_d)
  df_init$rmmd = df_init$mmd/ rowMeans(x_a)
  df_init$rldm = df_init$ldm/ rowMeans(x_a)
  df_init$rxbar_dm = df_init$xbar_dm/ rowMeans(x_a)
  df_init$rci_lower_z = df_init$ci_lower_z/ rowMeans(x_a)
  df_init$rci_upper_z = df_init$ci_upper_z/ rowMeans(x_a)
  df_init$rmu_dm = df_init$mu_dm/ rowMeans(x_a)
  
  df_list = list()
  # Sample mean and relative mean (errors computed for reference since statistic is uncontrolled)
  df_list[[1]] <- quant_error_rate(df_init = df_init, lower_name = "xbar_dm", upper_name = "xbar_dm" ,  
                              gt_name = "mu_dm", use_absolute = TRUE)
  df_list[[2]] <- quant_error_rate(df_init = df_init, lower_name = "rxbar_dm", upper_name = "rxbar_dm" ,  
                              gt_name = "rmu_dm", use_absolute = TRUE)
  
  # Confidence intervals of the mean, z distribution
  df_list[[3]] <- quant_error_rate(df_init = df_init, lower_name = "ci_lower_z", upper_name = "ci_upper_z" ,  
                          gt_name = "mu_dm", use_absolute = FALSE)
  df_list[[4]] <- quant_error_rate(df_init = df_init, lower_name = "rci_lower_z", upper_name = "rci_upper_z" ,  
                              gt_name = "mu_dm", use_absolute = FALSE)
  
  # Confidence intervals of the mean, z distribution
  df_list[[5]] <- quant_error_rate(df_init = df_init, lower_name = "ldm", upper_name = "mmd" ,  
                              gt_name = "mu_dm", use_absolute = TRUE)
  df_list[[6]] <- quant_error_rate(df_init = df_init, lower_name = "rldm", upper_name = "rmmd" ,  
                            gt_name = "rmu_dm", use_absolute = TRUE)
  df_err = do.call("cbind", df_list)

  return(df_err)
  
}
# 
# # Quantify coverage error of the mean (how often abs(x_bar) < abs(mu))
# abs_diff_xbar_mu               <- abs(rowMeans(x_d)) - abs(df$mu_d)
# df$mean_abs_diff_xbar_mu  <- mean(abs_diff_xbar_mu)
# n_errors                       <- sum(abs_diff_xbar_mu < 0)
# df$mean_xbar_error_rate   <- n_errors / n_samples
# df$pval_xbar_err_eq_zero  <- binom.test(
#   n_errors, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
# df$pval_xbar_err_eq_alpha <- binom.test(
#   n_errors, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
# # Difference rx_bar to rmu
# abs_diff_rxbar_mu               <- abs(rowMeans(x_d)/rowMeans(x_a)) -
#   abs(df$mu_d/df$mu_a)
# df$mean_abs_diff_rxbar_mu  <- mean(abs_diff_rxbar_mu)
# n_errors                        <- sum(abs_diff_rxbar_mu < 0)
# df$mean_rxbar_error_rate   <- n_errors / n_samples
# if (!is.nan(n_errors) && !is.na(n_errors) && !is.infinite(n_errors)) {
#   df$pval_rxbar_err_eq_zero  <- binom.test(
#     n_errors, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
#   df$pval_rxbar_err_eq_alpha <- binom.test(
#     n_errors, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
# } else {df$pval_rxbar_err_eq_zero= NaN; df$pval_rxbar_err_eq_alpha = NaN}
# 
# # Calculate the mmd from samples from the difference distribution
# mmd_d = apply(x_d, 1, function (x) mmd_normal_zdist(x, conf.level = 0.95) )
# # Difference mmd to mu
# #----------------------------------------------------------------------
# abs_diff_mmd_mu <- mmd_d - abs(df$mu_d)
# # Relative difference mmd to mu
# df$mean_diff_mmd_mu        <- mean(abs_diff_mmd_mu)
# df$mean_diff_mmd_mu_vmu    <- df$mean_diff_mmd_mu / abs(df$mu_d)
# df$mean_diff_mmd_mu_vsigma <- df$mean_diff_mmd_mu / df$sigma_d
# # Error rate mmd_d and mu
# n_errors                        <- sum(abs_diff_mmd_mu < 0)
# df$mean_mmd_error_rate     <- n_errors / n_samples
# # Caluate p-values for test against error rate == 0
# df$pval_mmd_err_eq_zero <- binom.test(
#   n_errors, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
# # Caluate p-values for test against error rate == alpha
# df$pval_mmd_err_eq_alpha <- binom.test(
#   n_errors, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
# 
# # Relative MMD: Difference r-mmd (mmd/ sample mean) to rmu
# #----------------------------------------------------------------------
# rmmd                         <- mmd_d / abs(rowMeans(x_a))
# df$mean_rmmd            <- mean(rmmd)
# # Quality check: no means of group a should be below zero for relative change
# df$fract_neg_x_bar_a    <- sum(rowMeans(x_a)<0) / n_samples
# diff_rmmd_rmu                <- rmmd - abs(df$rmu)
# df$mean_diff_rmmd_rmu   <- mean(diff_rmmd_rmu)
# # Error rate rmmd > rmu
# n_errors                     <- sum(diff_rmmd_rmu < 0)
# df$mean_rmmd_error_rate <- n_errors / n_samples
# # Calculate p-values for test rmmd error rate == 0
# if (!is.nan(n_errors) && !is.na(n_errors) && !is.infinite(n_errors)) {
#   df$pval_rmmd_err_eq_zero <- binom.test(
#     n_errors, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
#   # Calculate p-values for test rmmd error rate == alpha
#   df$pval_rmmd_err_eq_alpha <- binom.test(
#     n_errors, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
# }else {df$pval_rmmd_err_eq_zero= NaN; df$pval_rmmd_err_eq_alpha = NaN}
# 

quant_error_rate <-function(df_init, lower_name=NULL, upper_name=NULL, gt_name, 
                            error_rate_alpha = 0.05, binom.conf.level = 0.95,
                            use_absolute = TRUE){
  #' @description calculates coverage error for either single or dual bound 
  #' confidence intervals
  #' @param df_init input data frame where rows are samples and columns are 
  #' statistics and parameters
  #' @param lower_name string of column name of statistic used for lower bound of 
  #' parameter
  #' @param upper_name string of column name of statistic used for upper bound of 
  #' parameter
  #' @param gt_name string of parameter used as ground truth to compare bounds to
  #' @param error_rate_alpha significance level used to for statistics
  #' @param binom.conf.level confidence level for binomial tests
  #' @param use_absolute flag to compare magnitude of bounds and parameter
  #' @return dataframe of error rates, with column names tailored to the name of 
  #' the statistic and parameter used for the analysis

  # Number of samples (trials) for calculating candidate statistics, total for error rate
  n_samples = dim(df_init)[1]
   # Initialize 
  df = data.frame(0); colnames(df) <- paste("mean_abs_diff_", lower_name,"_to_", gt_name, sep="")
  if_binom_false <- function(x) {if(!x) {xp=2.2e-16} else {xp=x}}
  
  if (use_absolute) {abs_str = "abs_"}else{abs_str=""}
    
  # save(list = ls(all.names = TRUE), file = "temp/quant_error_rate.RData",envir = environment())
  # load(file = "temp/quant_error_rate.RData")
  
  
  # Error rate for lower bounds of parameter (how often lower bounds more than parameter)
  if (!is.null(lower_name)){
    
    # Compare magnitude of stat to parameter, or not
    if (use_absolute) {
    diffs_lower <- abs(df_init[[gt_name]]) - abs(df_init[[lower_name]])
    } else {diffs_lower <- df_init[[gt_name]] - df_init[[lower_name]]}
    # Calculate difference from stat to parameter (ground truth), and calculate error rate
    df[[paste("mean_", abs_str,"diff_", lower_name,"_to_", gt_name, sep="")]] <- mean(unname(diffs_lower))
    n_errors_lower                       <- sum(diffs_lower>0)
    df[[paste("mean_err_", abs_str, lower_name, "_gt_", gt_name, sep="")]] <- n_errors_lower / n_samples
    # Test if error rate is equal to zero or alpha
    if (is.finite(n_errors_lower)) {
      df[[paste("pval_err_eq_zero_",abs_str, lower_name, "_gt_", gt_name, sep="")]] <- 
        if_binom_false(binom.test(n_errors_lower, n_samples, p = 0.00, 
                                  alternative = "two.sided", conf.level = binom.conf.level)$p.value)
      df[[paste("pval_err_eq_alpha_",abs_str, lower_name, "_gt_", gt_name, sep="")]] <- 
        if_binom_false(binom.test(n_errors_lower, n_samples, p = error_rate_alpha, 
                                  alternative = "two.sided", conf.level = binom.conf.level)$p.value)
    } else {
      df[[paste("pval_err_eq_zero_",abs_str, lower_name, "_gt_", gt_name, sep="")]] <- NaN; 
      df[[paste("pval_err_eq_alpha_",abs_str, lower_name, "_gt_", gt_name, sep="")]] <- NaN
    }
  }
  
  # Error rate for upper bounds of parameter (how often upper bounds less than parameter)
  if (!is.null(upper_name)){
    # Calculate difference from stat to parameter (ground truth), and calculate error rate
    if (use_absolute) {
      diffs_upper <-   abs(df_init[[upper_name]]) - abs(df_init[[gt_name]]) 
    } else {diffs_upper <- df_init[[upper_name]]} -     df_init[[gt_name]]
    
    df[[paste("mean_", abs_str,"diff_", upper_name,"_to_", gt_name, sep="")]] <- mean(unname(diffs_upper))
    n_errors_upper                       <- sum(diffs_upper<0)
    df[[paste("mean_err_", abs_str, upper_name, "_lt_", gt_name, sep="")]] <- n_errors_upper / n_samples
    # Test if error rate is equal to zero or alpha
    if (is.finite(n_errors_upper)) {
      df[[paste("pval_err_eq_zero_",abs_str, upper_name, "_lt_", gt_name, sep="")]] <- 
        if_binom_false(binom.test(n_errors_upper, n_samples, p = 0.00, 
                                  alternative = "two.sided", conf.level = binom.conf.level)$p.value)
      df[[paste("pval_err_eq_alpha_",abs_str, upper_name, "_lt_", gt_name, sep="")]] <- 
        if_binom_false(binom.test(n_errors_upper, n_samples, p = error_rate_alpha,
                                  alternative = "two.sided", conf.level = binom.conf.level)$p.value)
    } else {
      df[[paste("pval_err_eq_zero_",abs_str, upper_name, "_lt_", gt_name, sep="")]] <- NaN; 
      df[[paste("pval_err_eq_alpha_",abs_str, upper_name, "_lt_", gt_name, sep="")]] <- NaN
    }
  }
  
  
  # Error rate for lower bounds and upper bounds of parameter (how often 
  # parameter is out of bounds)
  if (!is.null(lower_name) & !is.null(upper_name)){
    n_errors_both <- 
      sum((df_init[[lower_name]] > df_init[[gt_name]]) | 
            df_init[[upper_name]] < df_init[[gt_name]])
    df[[paste("mean_bound_err_", lower_name, "_", upper_name, sep="")]] <- n_errors_both / n_samples
    if (is.finite(n_errors_both)) {
      df[[paste("pval_bound_err_eq_zero_", upper_name, "_lt_", gt_name, sep="")]] <- 
        if_binom_false(binom.test(
          n_errors_both, n_samples, p = 0.00, alternative = "two.sided", 
          conf.level = binom.conf.level)$p.value)
      df[[paste("pval_bound_err_eq_alpha_", upper_name, "_lt_", gt_name, sep="")]] <-
        if_binom_false(binom.test(
          n_errors_both, n_samples, p = error_rate_alpha, alternative = "two.sided", 
          conf.level = binom.conf.level)$p.value)
    } else {
      df[[paste("pval_bound_err_eq_zero_", upper_name, "_lt_", gt_name, sep="")]] <- NaN; 
      df[[paste("pval_bound_err_eq_alpha_", upper_name, "_lt_", gt_name, sep="")]] <- NaN
    }
  }
  
  
  return(df)
  
}





row_locate_binary_bounds <- function (xb){
  #' Given a logical matrix, locates the border between true and false with a simple algroithm
  #' TODO|: replace algorithm, this one is not effective
  
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
  # browser()
  return(row_max_inds)
  
}


locate_bidir_binary_thresh <- function(ind_var = "mmd", mus = NULL, sigmas, n_samples, n_obs,
                                       temp_path, mu_ov_sigmas = NULL, rand.seed=0,
                                       overwrite = TRUE,  mus_a = 0, sigmas_a = 0,
                                       is_parallel_proc = TRUE){
  #' @description Calculates the coverage error of x_bar, rx_bar, mmd, and rmmd
  #'
  #' @param ind_var statistic evaluted for coeverage error
  #' @param mus vector of range of population mu values to be quantified
  #' @param sigmas vector of range of population sigma values to be quantified
  #' @param n_samples number of samples drawn or each combination of mu and sigma
  #' @param n_obs number of observations drawn for each simulation
  #' @param temp_path path to store temp data
  #' @param mu_ov_sigmas vector of range of mu/sigma values to be quantified, 
  #' @param rand.seed seed for random number generation
  #' @param overwrite flag to overwrite temp data (if not, load temp if exist)
  #' @param mus_a vector of means for the control group (for quantifying coverage
  #'  error of relative x_bar and mmd, the mean of of the control group must be
  #'  specified.
  #' @param sigmas_a vector of stds for the control group (for quantifying coverage
  #'  error of relative x_bar and mmd, the std of the control group must be
  #'  specified.
  #' @param is_parallel_proc flag whether simulations are run in parallel across cores
  #' 
  #' @return df_crit dataframe that identifies the boundary between error regions

  # Create temp directory and set random seet
  dir.create(temp_path, showWarnings = FALSE, recursive = TRUE)
  set.seed(rand.seed)
  
  save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  # load(file = "temp/debug.RData")
  # Create temp dir if it does not exist
  
  
  #' Locate all 4 error transition boundaries
  #' Positive and negative direction, 0 and alpha error
  
  n_cols <- max(c(length(mus), length(mu_ov_sigmas)))
  if (!is.null(mus)) {
    col_values = mus; 
    varname <- "critical_mu"; 
    neg_mu_ov_sigmas = mu_ov_sigmas;
  } else {
    col_values = mu_ov_sigmas; 
    varname <- "critical_mu_over_sigma"; 
    neg_mu_ov_sigmas=-mu_ov_sigmas

  }
  

  # Run simulations calculating error of mmd with mu and sigma swept
  df_right <- quant_coverage_errors(mus_ao = mus, sigmas_ao = sigmas, 
                                n_samples = n_samples, n_obs = n_obs, paste(temp_path,"_right.rds"),
                                mu_vsigmas_ao = mu_ov_sigmas, overwrite = overwrite,
                                mus_a = mus_a, sigmas_a = sigmas_a, is_parallel_proc = is_parallel_proc) 
  df_right$side <- as.factor("Right")
  
  # Run simulations calculating error of mmd with mu and sigma swept
  df_left <- quant_coverage_errors(mus_ao = mus, sigmas_ao = sigmas, 
                                      n_samples = n_samples, n_obs = n_obs, 
                                   paste(temp_path,"_right.rds"),
                               mu_vsigmas_ao= neg_mu_ov_sigmas, overwrite = overwrite,
                               mus_a = mus_a, sigmas_a = sigmas_a, is_parallel_proc = is_parallel_proc) 
  df_left$side <- as.factor("Left")
  

  ind_var_zero <- paste("pval_",ind_var,"_err_eq_zero",sep="")
  ind_var_alpha <- paste("pval_",ind_var,"_err_eq_alpha",sep="")
  
  # Equivalence test versus middle column of same row
  p_threshold = 0.05 /(n_cols*length(sigmas))
  df_zero <- rbind(
    tibble(er = "0", side = "right", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_right[[ind_var_zero]]  < p_threshold)),
    tibble(er = "alpha", side = "right", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_right[[ind_var_alpha]]   < p_threshold)))
  # Convert index of error rate transition to mu/sigma value
  df_zero[[varname]] <- 
    approx(x=1:n_cols,y = abs(col_values),
           xout = df_zero$critical_val_ind, 
           n = length(n_cols)*2L-1L)$y
  
  df_alpha <- rbind(
    tibble(er = "0", side = "left", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_left[[ind_var_zero]]    < p_threshold)),
    tibble(er = "alpha", side = "left", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_left[[ind_var_alpha]]   < p_threshold)))
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



