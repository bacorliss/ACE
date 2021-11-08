
# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(foreach)
p_load(doParallel)
p_load(tidyr)
p_load(stringr)
source("R/row_stats_toolbox.R")
source("R/parallel_utils.R")
source("R/aces.R")

mdm_functions <- parse_functions_source("R/aces.R")
# ldm_functions <- parse_functions_source("R/ldm.R")
row_effect_size_functions <- parse_functions_source("R/row_stats_toolbox.R")
# rationormal_functions <- parse_functions_source("R/rationormal_toolbox.R")

RowVar <- function(x, ...) rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
# Test if a value is within an interval (inclusive)
within_ci <- function(ci, x) x >= ci[1] & x <= ci[2]


quant_coverage_errors <-
  function(mus_a, sigmas_a, n_a, mus_b, sigmas_b, n_b, alphas, n_samples, out_path, 
           mu_vsigmas_dm = NA, overwrite = FALSE, rand.seed = 0,
           is_parallel_proc = TRUE, raw_error = TRUE, rel_error = TRUE, included_stats = NULL) {
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
    #'  @param mu_vsigmas_dm used in place of the mus_ao vector when looking at 
    #'  normalized values of mu
    #'  
    #'  @return df dataframe that stores population parameters and statistics
    
    save(list = ls(all.names = TRUE), file = "temp/quant_coverage_errors.RData",envir = environment())
    # load(file = "temp/quant_coverage_errors.RData")
    dir.create(dirname(out_path),showWarnings = FALSE, recursive = TRUE)
    
    # For any pop param not equal in length to n_sims, expand
    n_mus = max(length(mus_a), length(mus_b), length(mu_vsigmas_dm))
    n_sigmas = max(length(sigmas_a), length(sigmas_b))
    if (length(mus_a)==1) {mus_a = rep(mus_a, n_mus)}
    if (length(mus_b)==1) {mus_b = rep(mus_b, n_mus)} 
    if (length(sigmas_a)==1) {sigmas_a = rep(sigmas_a, n_sigmas)}
    if (length(sigmas_b)==1) {sigmas_b = rep(sigmas_b, n_sigmas)}
    if (length(mu_vsigmas_dm)==1) {mu_vsigmas_dm = rep(mu_vsigmas_dm, n_mus)}
    
    mus_dm = mus_b - mus_a
    sigmas_dm = sqrt(sigmas_a^2/n_a + sigmas_b^2/n_b) 
    
    # Store results to disk since calculations are significant
    set.seed(rand.seed)
    if (!file.exists(out_path) || overwrite) {
      
      # Matrix diff and error of mdm
      if (all(is.na(mu_vsigmas_dm))) {dimnames = list(sigmas_dm,mus_dm)
      } else {dimnames = list(sigmas_dm, mu_vsigmas_dm)}
      
      
      # Initialize matrix in data form
      param_col_list <- c("mu_a", "sigma_a","n_a", "mu_b","sigma_b", "n_b", 
                          "alpha", "mu_dm", "sigma_dm", "mu_vsigma_dm",
                          "rmu_dm", "n_samples", "alpha")
      
      init_mat <- matrix(NA, nrow = n_sigmas, ncol = n_mus, dimnames = dimnames)
      df_mat <- tibble(mu_a = init_mat);
      for (n in seq_along(param_col_list)) {df_mat[[param_col_list[n]]] <- init_mat }
      
      for (r in seq(1, n_sigmas, by = 1)) {
        # Calculate mu if only mu/sigma and sigma are known
        if (!all(is.na(mu_vsigmas_dm))) {
          mus_dm = mu_vsigmas_dm * sigmas_dm[r]
          mus_b = mus_a + mus_dm
        }
        
        for (c in seq(1, n_mus, by = 1)) {
          df_mat$mu_a[r,c]     <- mus_a[c]
          df_mat$sigma_a[r,c]  <- sigmas_a[r]
          df_mat$n_a[r,c]  <- n_a
          
          df_mat$mu_b[r,c]     <- mus_b[c]
          df_mat$sigma_b[r,c]  <- sigmas_b[r]
          df_mat$n_b[r,c]  <- n_b
          
          # Difference in means group
          df_mat$mu_dm[r,c]     <- mus_dm[c]
          df_mat$sigma_dm[r,c]  <- sigmas_dm[r]
          df_mat$alpha[r,c]  <- alphas
          
          df_mat$mu_vsigma_dm[r,c] <- df_mat$mu_dm[r,c]/df_mat$sigma_dm[r,c]
          df_mat$rmu_dm[r,c] <- df_mat$mu_dm[r,c] / df_mat$mu_a[r,c]
          df_mat$n_samples[r,c] = n_samples
          
        }
      } 
      
      # Linearize matrix dataframe for parallel processing
      df_lin <- tibble(mu_a = as.vector(init_mat));
      for (n in seq_along(param_col_list)) { df_lin[[param_col_list[n]]] <- 
        as.vector(df_mat[[param_col_list[n]]]) }
      
      
      # Process parallel or serially
      if (is_parallel_proc) { 
        print("Starting parallel cluster...")
        # browser()
        cl = makeCluster(max(c(floor(detectCores()[1]*.9), 1)))
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
        
        df_lin2 <- 
          foreach(n = seq(1,n_sigmas*n_mus,1),
                  .export = c(mdm_functions, row_effect_size_functions,
                              "quant_coverage_error","quant_error_rate"), 
                  .combine = rbind, .packages = c("tidyr", "cubature")) %dopar% {
                    #calling a function
                    tempMatrix <- quant_coverage_error(df = df_lin[n,], 
                                                       raw_error = raw_error,
                                                       rel_error = rel_error,
                                                       included_stats = included_stats) 
                    tempMatrix
                  }
        # stopCluster(cl)
      } else {            # Process effect sizes serially
        # browser();
        df_list <- list()
        for (n in seq(1,n_sigmas*n_mus,1)) {
          print(n)
          df_list[[n]] <- quant_coverage_error(df_lin[n,], 
                                               raw_error = raw_error,
                                               rel_error = rel_error,
                                               included_stats = included_stats) 
          # browser();
        }
        df_lin2 <- do.call("rbind", df_list)
      }
      
      # Convert linear data frame to dataframe of matrices
      df_mat2 <- df_mat; col_list <- colnames(df_lin2)
      for (n in seq_along(col_list)) { 
        df_mat2[[col_list[n]]] <- matrix(df_lin2[[col_list[n]]], nrow = n_sigmas, 
                                         ncol = n_mus, dimnames = dimnames)
      }

      
      # Save an object to a file
      saveRDS(df_mat2, file = out_path)
    } else {
      # Restore the object
      df_mat2 <- readRDS(file = out_path)
    }
    # browser();
    save(list = ls(all.names = TRUE), file = "temp/quant_coverage_errors.RData",
         envir = environment())
    return(df_mat2)
  }



quant_coverage_error <-  function(df, raw_error = TRUE, rel_error = TRUE, verbose = FALSE, 
                                  rand.seed = NULL, enforce_grand_mean = FALSE, 
                                  included_stats = NULL) {
  #' @description Calculates the coverage error of x_expar, rx_expar, mdm, and rmdm
  #'
  #' @param df: a single row dataframe generated from agreement_contest
  #'  library, returned from generate_population_configs()
  #' @param n_samples: number of samples drawn to evaluate converage error
  #' @param n_obs: 
  #' @param enforce_grand_mean flag to force the grand mean (mean across all 
  #' samples and observations) to be exactly equal to the population mean. This 
  #' is used for debugging/ research purposes and not for actual simulations.
  #' @return null, exports figures to disk
  
  save(list = ls(all.names = TRUE), file = "temp/quant_coverage_error.RData",envir = environment())
  # load(file = "temp/quant_coverage_error.RData")
  
  if (is.null(included_stats)) {stop("quant_coverage_error: arg included_stats cannot be null")}
  if (verbose) {print(sprintf("A: %f, B: %f", df$mu_a, df$mu_b))}
  
  if (!is.null(rand.seed)) {set.seed(rand.seed)}
  
  # Control group (For use with two sample cases)
  x_ctrl <- matrix(rnorm(df$n_samples * df$n_a, df$mu_a, df$sigma_a), ncol = df$n_a)
  # Difference group (for simplicity experiment sample not calculated)
  x_exp <- matrix(rnorm(df$n_samples * df$n_b, df$mu_b, df$sigma_b), ncol = df$n_b)
  # Each row is a separate sample, columns are observations
  
  if (enforce_grand_mean) {
    # For debugging purposes, forces grandmean of all samples to population mu
    x_ctrl <- x_ctrl - (mean(x_ctrl) - df$mu_a)
    x_exp <-  x_exp -  (mean(x_exp) -  df$mu_b)
  }
  
  # Only compute requested statistics
  df_include = data.frame(matrix(rep(TRUE,length(included_stats)), nrow = 1, 
                                 dimnames = list(NULL, included_stats)))
  # browser()
  # Initial sample metrics
  # ci_mean = row_ci_mean_2s_zdist(m_c = x_ctrl, m_e = x_exp)
  
  # browser();
  
  df_init = data.frame(xbar_dm = rowMeans(x_exp) - rowMeans(x_ctrl))
  
  # Initialize list of comparison error each element has results with different 
  # metrics and groundtruths  
  df_list = list()
  
  # Raw error groundtruth
  df_init$mu_dm = df$mu_dm
  
  if (!is.null(df_include$xbar_dm)) { # raw difference in sample means
  df_list[[1]] <- quant_error_rate(df_init = df_init, lower_name = "xbar_dm", upper_name = "xbar_dm" ,  
                                   gt_name = "mu_dm", use_absolute = TRUE)
  }
  
  if (!is.null(df_include$ci_dm)) { # raw confidence interval of difference in means
  # Confidence intervals of the mean, z distribution
  df_list[[2]] <- quant_error_rate(df_init = df_init, lower_name = "ci_lower_z", upper_name = "ci_upper_z" ,
                                   gt_name = "mu_dm", use_absolute = FALSE)
  }
  
  
  
  # browser();
  
  if (!is.null(df_include$mdm)) {
    # aquantile_mdm = quantile(abs(rowMeans(x_exp) - rowMeans(x_ctrl)), 1-df$alpha[1])
    
    df_init$mdm = row_mdm(m_c = x_ctrl, m_e = x_exp, conf.level = 1 - df$alpha)
    df_list[[length(df_list)+1]] <- quant_error_rate(df_init = df_init, lower_name = NULL, upper_name = "mdm",
                                                     gt_name = "mu_dm", use_absolute = TRUE)
    # df_list[[length(df_list)]]$diff_mdm_aquantile_mdm = mean(df_init$mdm) - mean(df_init$aquantile_mdm)
  }
  
  if (!is.null(df_include$ldm)) {
    df_init$ldm = row_ldm_2s_zdist(x_ctrl, x_exp, conf.level = 1 - df$alpha)
    df_list[[length(df_list)+1]] <- quant_error_rate(df_init = df_init, lower_name = "ldm", upper_name = NULL,
                                                     gt_name = "mu_dm", use_absolute = TRUE)
  }
  
  
  if (!is.null(df_include$macb)) { 
    df_init$macb = row_macb_tdist_2sample(m_c = x_ctrl, m_e = x_exp, conf.level = 1 - df$alpha)
    df_list[[length(df_list)+1]] <- quant_error_rate(df_init = df_init, lower_name = NULL, upper_name = "macb",
                                                     gt_name = "mu_dm", use_absolute = TRUE)
  }
  
  
  if (!is.null(df_include$macb_aov2)) {   
    df_init$macb_aov2 = row_macb_tdist_2sample(m_c = x_ctrl, m_e = x_exp, conf.level = 1 - 0.5*df$alpha)
    df_list[[length(df_list)+1]] <- quant_error_rate(df_init = df_init, lower_name = NULL, upper_name = "macb_aov2",
                                                     gt_name = "mu_dm", use_absolute = TRUE)
  }
  
  
  # Ground truth for relative difference in means (population values)
  df_init$rmu_dm = (df$mu_b - df$mu_a)/df$mu_a
  # Relative difference in sample means
  df_init$rxbar_dm = (rowMeans(x_exp) - rowMeans(x_ctrl))/ rowMeans(x_ctrl)
  
  if (!is.null(df_include$rxbar_dm)) {
  df_list[[length(df_list)+1]] <-
    quant_error_rate(df_init = df_init, lower_name = "rxbar_dm", upper_name = "rxbar_dm" ,
                     gt_name = "rmu_dm", use_absolute = TRUE)
  }
  
  if (!is.null(df_include$rmdm)) {
    # aquantile_rmdm = quantile(abs(rowMeans(x_exp) - rowMeans(x_ctrl))/abs(rowMeans(x_ctrl)), 1-df$alpha[1])
    df_init$rmdm = row_rmdm(x_ctrl, x_exp, conf.level = 1-df$alpha, mdms = df_init$mdm)
    if (any(is.nan(df_init$rmdm))) { browser() }
    df_list[[length(df_list) + 1]] <-
      quant_error_rate(df_init = df_init, lower_name = NULL, upper_name = "rmdm" ,
                       gt_name = "rmu_dm", use_absolute = TRUE)
    
  }
  
  if (!is.null(df_include$rldm)) {
    df_list[[length(df_list) + 1]] <-
      quant_error_rate(df_init = df_init, lower_name = "ldm", upper_name = NULL,
                       gt_name = "rmu_dm", use_absolute = TRUE)
  }
  
  if (!is.null(df_include$rci_dm)) {
    df_init$rci_lower_z = df_init$ci_lower_z/ rowMeans(x_ctrl)
    df_init$rci_upper_z = df_init$ci_upper_z/ rowMeans(x_ctrl)
    df_list[[length(df_list) + 2]] <-
      quant_error_rate(df_init = df_init, lower_name = "rci_lower_z", upper_name = "rci_upper_z" ,
                       gt_name = "mu_dm", use_absolute = FALSE)
  }
  
  
  if (!is.null(df_include$rxbar_eoc)) {
    df_init$rxbar_eoc = row_ratio_normal_eoc(x_ctrl, x_exp, conf.level = 1 - df$alpha)
    df_init$rmu_eoc = df$mu_b/df$mu_a
    df_list[[length(df_list) + 1]] <-
      quant_error_rate(df_init = df_init, lower_name = NULL, upper_name = "rxbar_eoc" ,
                       gt_name = "rmu_eoc", use_absolute = TRUE)
  }
  
  
  
  # browser();
  df_err = do.call("cbind", df_list)
  return(df_err)
  
  
}

  
  


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
    
    df[[paste("mean_", abs_str, lower_name, sep="")]] = mean(df_init[[lower_name]])
    
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
    
    df[[paste("mean_", abs_str, upper_name, sep="")]] = mean(df_init[[upper_name]])
    
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





locate_bidir_binary_thresh <- 
  function(ind_var = "mdm", pop_var = "mu_dm", mus_a, sigmas_a, n_a, 
           mus_b, sigmas_b, n_b, mu_vsigmas_dm = NA, alphas = 0.05, n_samples, 
           temp_path, overwrite = overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE)
  {
    #' @description Calculates the coverage error of x_expar, rx_expar, mdm, and rmdm
    #'
    #' @param ind_var statistic evaluted for coeverage error
    #' @param mus vector of range of population mu values to be quantified
    #' @param sigmas vector of range of population sigma values to be quantified
    #' @param n_samples number of samples drawn or each combination of mu and sigma
    #' @param n_obs number of observations drawn for each simulation
    #' @param temp_path path to store temp data
    #' @param mu_vsigmas_dm vector of range of mu/sigma values to be quantified, 
    #' @param rand.seed seed for random number generation
    #' @param overwrite flag to overwrite temp data (if not, load temp if exist)
    #' @param mus_a vector of means for the control group (for quantifying coverage
    #'  error of relative x_expar and mdm, the mean of of the control group must be
    #'  specified.
    #' @param sigmas_a vector of stds for the control group (for quantifying coverage
    #'  error of relative x_expar and mdm, the std of the control group must be
    #'  specified.
    #' @param is_parallel_proc flag whether simulations are run in parallel across cores
    #' 
    #' @return df_crit dataframe that identifies the boundary between error regions
    
    # Create temp directory and set random seet
    dir.create(temp_path, showWarnings = FALSE, recursive = TRUE)
    set.seed(rand.seed)
    
    # save(list = ls(all.names = TRUE), file = "temp/locate_bidir_binary_thresh.RData",envir = environment())
    # load(file = "temp/locate_bidir_binary_thresh.RData")
    # Create temp dir if it does not exist
    
    # Difference in means
    mus_dm = mus_b - mus_a
    sigmas_dm = sqrt(sigmas_a^2/n_a + sigmas_b^2/n_b) 
    
    
    #' Locate all 4 error transition boundaries
    #' Positive and negative direction, 0 and alpha error
    
    n_cols <- max(c(length(mus_dm), length(mu_vsigmas_dm)))
    if (!any(is.na(mus_dm)) & !any(is.null(mus_dm))) {
      col_values = mus_dm; 
      varname <- "critical_mu"; 
      neg_mu_vsigmas_dm = -mu_vsigmas_dm;
    } else {
      col_values = mu_vsigmas_dm; 
      varname <- "critical_mu_over_sigma"; 
      neg_mu_vsigmas_dm=-mu_vsigmas_dm
      
    }
    
    
    # browser();
    
    # Run simulations calculating error of mdm with mu and sigma swept
    df_right <- quant_coverage_errors( mus_a, sigmas_a, n_a, mus_b, sigmas_b, n_b, alphas, n_samples, 
                                       str_replace(temp_path, ".rds$","_right.rds"), 
                                       mu_vsigmas_dm = mu_vsigmas_dm, overwrite = overwrite, rand.seed=0,
                                       is_parallel_proc = is_parallel_proc, 
                                       raw_error = raw_error, rel_error = rel_error,
                                       included_stats = c(ind_var)) 
    df_right$side <- as.factor("Right")
    
    # Run simulations calculating error of mdm with mu and sigma swept
    df_left <- quant_coverage_errors( mus_a, sigmas_a, n_a, mus_b, sigmas_b, n_b, alphas, n_samples, 
                                      str_replace(temp_path, ".rds$","_right.rds"), 
                                      mu_vsigmas_dm = neg_mu_vsigmas_dm, overwrite = overwrite, rand.seed=0,
                                      is_parallel_proc = is_parallel_proc, 
                                      raw_error = raw_error, rel_error = rel_error,
                                      included_stats = c(ind_var)) 
    df_left$side <- as.factor("Left")
    
    save(list = ls(all.names = TRUE), file = "temp/locate_bidir_binary_thresh.RData",envir = environment())
    # load(file = "temp/locate_bidir_binary_thresh.RData")
    # browser();
    
    
    ind_var_zero <-  paste("pval_err_eq_zero_abs_",ind_var,"_lt_", pop_var,sep="")
    ind_var_alpha <- paste("pval_err_eq_alpha_abs_",ind_var,"_lt_", pop_var,sep="")
    
    # Equivalence test versus middle column of same row
    p_threshold = 0.05 /(n_cols*length(sigmas_dm))
    df_right_border <- rbind(
      tibble(er = "0", side = "right", sigma = sigmas_dm, 
             critical_val_ind = row_locate_binary_bounds(
               df_right[[ind_var_zero]]  < p_threshold,  true_side = "right")),
      tibble(er = "alpha", side = "right", sigma = sigmas_dm, 
             critical_val_ind = row_locate_binary_bounds(
               df_right[[ind_var_alpha]]   < p_threshold, true_side = "left")))
    # Convert index of error rate transition to mu/sigma value
    df_right_border[[varname]] <- 
      approx(x=1:n_cols,y = abs(col_values),
             xout = df_right_border$critical_val_ind, 
             n = length(n_cols)*2L-1L, rule = 2)$y
    
    df_left_border <- rbind(
      tibble(er = "0", side = "left", sigma = sigmas_dm, 
             critical_val_ind = row_locate_binary_bounds(
               df_left[[ind_var_zero]]    < p_threshold, true_side = "right")),
      tibble(er = "alpha", side = "left", sigma = sigmas_dm, 
             critical_val_ind = row_locate_binary_bounds(
               df_left[[ind_var_alpha]]   < p_threshold, true_side = "left")))
    # Convert index of error rate transition to mu/sigma value
    df_left_border[[varname]] <- 
      approx(x=1:n_cols,y = abs(col_values),
             xout = df_left_border$critical_val_ind, 
             n = length(n_cols)*2L-1L, rule = 2)$y
    
    # Concatenate right and left dataframes of results
    df_crit <- rbind(df_right_border,df_left_border)
    
    df_crit[[varname]][df_crit$side=="left"] = -
      abs(df_crit[[varname]][df_crit$side=="left"])
    # * abs() added just in case we are running the code multiple times when debugging
    
    # browser()
    
    return(df_crit)
    
  }





row_locate_binary_bounds <- function (xa, true_side = "left") {
  #' Given a logical matrix, locates the border between true and false with a simple algroithm
  #' TODO|: replace algorithm, this one is not effective
  

  save(list = ls(all.names = TRUE), file = "temp/row_locate_binary_bounds.RData",envir = environment())
  # load(file = "temp/row_locate_binary_bounds.RData")
  # browser()
  
  # Assumptions location of TRUE and false regions known and oriented correctly
  # TRUE: from left to right
  # FALSE: from right to left
  # Boundary location is between the point index (0.5, 1.5, 2.5 etc.)
  
  
  if ((true_side != "left") & (true_side != "right")) 
    {stop('row_locate_binary_bounds: uknown orientation option, see help')}
  
  if (true_side == "right") {xb = xa[,seq(dim(xa)[2],1,-1)]} else {xb = xa}
  
  row_max_inds = rep(0, dim(xb)[1])
  
  for (i in seq(1, dim(xb)[1])) {
    
    x = unname(xb[i,])
    # Make vectors of TP, FP, TN, FN
    TP =  c(0,cumsum(x))
    FP = c(rev(cumsum(rev(x))),0)
    
    TN = c(rev(cumsum(rev(!x))),0)
    FN = c(0, cumsum(!x))
    # rbind(TP, FP, TN, FN)
    
    TPR = TP / (TP + FN) # Sensitivity, Recall
    FPR = FP / (FP + TN)
    TNR = TN / (FP + TN) # Specificity
    
    PREC = TP / (TP + FP)
    REC = TP / (TP + FN)
    F1_SCORE = 2 * (PREC * REC) / (PREC + REC)
    
    # Adjust index to between datapoints in F1_score
    thresh = which.max(F1_SCORE) - .5
    
    # 
    if (length(thresh)!=0) {row_max_inds[i] = thresh
    } else if (all(!x)) { row_max_inds[i] = 0         +  0.5 
    } else if (all( x)) { row_max_inds[i] = length(x) +  0.5 
    }
    
  }
  
  if (true_side == "right") {
    row_max_inds = length(x) - row_max_inds + 1 
  }
    
  return(row_max_inds)
  
}



error_test_codes <-function(is_error_rate_zero, is_error_rate_alpha) {
  #' @description Assign codes for error rate whether null hypothesis is rejected
  #' for use in a heatmap of hypothesis outcomes. Test if proportion of mdms that
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

