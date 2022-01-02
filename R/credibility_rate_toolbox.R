# Quantifies coverage error of the mdm (defined as how often abs(mu)>mdm_95 
# from repeated samples)
# Results computed in a grid with mu and sigma swept on each exis respectively

# Load required packages
#______________________________________________________________________________
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# Load packages
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(broom)
p_load(gridExtra)
p_load(grid)
p_load(rlang)
p_load(colorspace)
p_load(VGAM)
p_load(boot)
p_load(dplyr)
p_load(cowplot)
p_load(tidyverse)
p_load(VGAM)
p_load(gplots)
p_load(RColorBrewer)
p_load(tidyr)
p_load(docstring)
p_load(foreach)
p_load(doParallel)
# User defined functions
source("R/aces.R")
source("R/row_stats_toolbox.R")
# source("R/coverage_error_toolbox.R")


# Figure parameters
#______________________________________________________________________________
base_dir = "mdm_t"
# Script Parameters
fig_num = "4"
fig_path = file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE,recursive = TRUE)
n_samples <- 1e3
rand.seed <- 0
overwrite <- TRUE
is_parallel_proc <- TRUE







quant_cred_interval <-  function(df_sample, stat_name, n_mus = 25, n_sigmas = 27,
                                 p_thresh_mu = 1e-8, p_thresh_sigma = 1e-8, 
                                 method = "montecarlo")  {
  #' @description Calculates the coverage error of x_expar, rx_expar, mdm, and rmdm
  #'
  #' @param df: a single row dataframe generated from agreement_contest
  #'  library, returned from generate_population_configs()
  #' @param n_samples: number of samples drawn to evaluate converage error
  #' @param n_obs: 
  #' @return null, exports figures to disk
  # Quick run
  # n_mus = 100; n_sigmas = 100; p_thresh_mu = 1e-6; p_thresh_sigma = 1e-6
  save(list = ls(all.names = TRUE), file = "temp/quant_cred_interval.RData",envir = environment())
  # load(file = "temp/quant_cred_interval.RData")
  
 
  xbar <- -df_sample$xbar_a; ybar <- -df_sample$xbar_b
  m <- df_sample$n_a; n <- df_sample$n_b
  s2x <- df_sample$sd_a^2; s2y <- df_sample$sd_b^2
  shape1 <- .5*(m-1)
  scale1 <- .5*(m-1)*s2x
  shape2 <- .5*(n-1)
  scale2 <- .5*(n-1)*s2y
  
  
  if (method == "montecarlo") { 

    # Simulate mu and sigma values based on sample mean and std
    ss1Sims <- 1/rgamma(n = df_sample$n_samples, shape = shape1, rate = scale1)
    ss2Sims <- 1/rgamma(n = df_sample$n_samples, shape = shape2, rate = scale2)
    mu1Sims <- rnorm(n = df_sample$n_samples, mean = xbar, sd = sqrt(ss1Sims/m))
    mu2Sims <- rnorm(n = df_sample$n_samples, mean = ybar, sd = sqrt(ss2Sims/n))
    
    # Calculate percentile of mdm or rmdm of quantity
    if (stat_name=="mdm") {
      x <- mu1Sims - mu2Sims
      df_out <- tibble(cred_rate = sum(abs(x) < df_sample$mdm)/df_sample$n_samples, 
                       n_samples = df_sample$n_samples)
    } else if (stat_name=="rmdm") {
      x <- (mu2Sims - mu1Sims)/mu1Sims
      df_out <- tibble(cred_rate = sum(abs(x) < df_sample$rmdm)/df_sample$n_samples, 
                       n_samples = df_sample$n_samples)
    } else {stop("quant_cred_interval:unknown stat_name")}
    
  }else if (method=="discrete_pdf") {
    

    # Calculate range for sigma based on a threshold for pdf
    # Since distribution is for sigma2 1 and 2 separately, just do sigma2_1 
    # since we are making variance equal
    sigma2_1_lo.bound <- qgamma( p_thresh_sigma,     shape = shape1, rate = scale1)
    sigma2_1_hi.bound <- qgamma( 1 - p_thresh_sigma, shape = shape1, rate = scale1)
    sigma2_1_diff <- (sigma2_1_hi.bound - sigma2_1_lo.bound)/(n_sigmas-1)
    sigmas2_1 <- seq(sigma2_1_lo.bound, sigma2_1_hi.bound, sigma2_1_diff)
    # Probability that each bin of sigma will be drawn
    p_sigmas <- pgamma(sigmas2_1+sigma2_1_diff/2, shape = shape1, rate = scale1) - 
              pgamma(sigmas2_1-sigma2_1_diff/2, shape = shape1, rate = scale1)
    sigmas2_dm <- sigmas2_1/m + sigmas2_1/n
    
    # Calculate mu range for each sigma
    mu_lo.bound <- qnorm(p_thresh_mu,mean = df_sample$xbar_dm, sd = sqrt(sigmas2_dm))
    mu_hi.bound <- qnorm(1-p_thresh_mu,mean = df_sample$xbar_dm, sd = sqrt(sigmas2_dm))
    mu_diffs <- (mu_hi.bound - mu_lo.bound)/(n_mus-1)
    # mus <-  seq(mu_lo.bound, mu_hi.bound, mu_diffs)
    mus_2d <- t(mapply(function(x,y,z) seq(x,y,z), mu_lo.bound, mu_hi.bound,mu_diffs))
    p_mus = list()
    for (n in seq(1, n_sigmas)) {
      p_mus[[n]] <- 
        pnorm(mus_2d[n,] + mu_diffs[n]/2, mean = df_sample$xbar_dm, sd = sqrt(sigmas2_dm[n])) -
        pnorm(mus_2d[n,] - mu_diffs[n]/2, mean = df_sample$xbar_dm, sd = sqrt(sigmas2_dm[n]))
      
    }
    p_mus_2d <- do.call(rbind,p_mus)
    

    p_tot_list  = list()
    for (n in seq(1, n_sigmas)) {
     p_tot_list[[n]] = p_sigmas[n] * p_mus_2d[n,]
    }
    p_tot_2d <- do.call(rbind,p_tot_list)
    
    cred_rate <- sum(p_tot_2d * (mus_2d <df_sample[[stat_name]]))
  
    
    df_out <- df_sample
    df_out[[paste("n_trials_", stat_name , "_pass",sep="")]] = cred_rate 
    df_out[[paste("n_trials_", stat_name , "_fail",sep="")]] = 1-cred_rate 
    df_out[[paste("credibility_rate_", stat_name ,"_lt_",stat_name ,sep="")]] = cred_rate 
    
  }
  
  return(df_out)
}




process_cred_intervals <-function(xbars_a, sds_a, n_a, xbars_b, sds_b, n_b, alphas = 0.05,
                                  n_samples = n_samples, out_path,
                                  overwrite = TRUE, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                                  stat_name, method) {
  save(list = ls(all.names = TRUE), file = "temp/process_cred_intervals.RData",envir = environment())
  # load(file = "temp/process_cred_intervals.RData")
  dir.create(dirname(out_path),showWarnings = FALSE, recursive = TRUE)
  
  # # Only compute requested statistics
  # df_include = data.frame(matrix(rep(TRUE,length(included_stats)), nrow = 1, 
  #                                dimnames = list(NULL, included_stats)))
  
  # For any pop param not equal in length to n_sims, expand
  n_xbars = max(length(xbars_b), length(xbars_a))
  n_sds = max(length(sds_a), length(sds_b))
  if (length(xbars_b)==1) {xbars_b = rep(xbars_b, n_xbars)}
  if (length(xbars_a)==1) {xbars_a = rep(xbars_a, n_xbars)} 
  if (length(sds_a)==1) {sds_a = rep(sds_a, n_sds)}
  if (length(sds_b)==1) {sds_b = rep(sds_b, n_sds)}
  
  # Calculate difference in mean statistics
  xbars_dm = xbars_b - xbars_a
  sds_dm = sqrt(sds_a^2/n_a + sds_b^2/n_b) 
  
  # Store results to disk since calculations are significant
  set.seed(rand.seed)
  if (!file.exists(out_path) || overwrite) {
    
    # Matrix diff and error of mdm
    dimnames = list(sds_dm,xbars_dm)
    # Initialize matrix in data form
    param_col_list <- c("xbar_a", "sd_a","n_a", "xbar_b","sd_b", "n_b", 
                        "xbar_dm", "sd_dm", "rxbar_dm", "n_samples" ,"alpha")
    
    init_mat <- matrix(NA, nrow = n_sds, ncol = n_xbars, dimnames = dimnames)
    df_mat <- tibble(xbar_a = init_mat);
    for (n in seq_along(param_col_list)) {df_mat[[param_col_list[n]]] <- init_mat }
    
    for (r in seq(1, n_sds, by = 1)) {
      for (c in seq(1, n_xbars, by = 1)) {
        df_mat$xbar_a[r,c]     <- xbars_a[c]
        df_mat$sd_a[r,c]  <- sds_a[r]
        df_mat$n_a[r,c]  <- n_a
        
        df_mat$xbar_b[r,c]     <- xbars_b[c]
        df_mat$sd_b[r,c]  <- sds_b[r]
        df_mat$n_b[r,c]  <- n_b
        
        # Difference in means group
        df_mat$xbar_dm[r,c]     <- xbars_b[c] - xbars_a[c]
        df_mat$sd_dm[r,c]  <- sqrt(df_mat$sd_a[r,c]^2/df_mat$n_a[r,c] + 
                                     df_mat$sd_b[r,c]^2/df_mat$n_b[r,c]) 
        df_mat$alpha[r,c]  <- alphas
        
        df_mat$rxbar_dm[r,c] <- df_mat$xbar_dm[r,c] / df_mat$xbar_a[r,c]
        df_mat$n_samples[r,c] = n_samples
        
        
      }
    } 
    
    
    
    # Linearize matrix dataframe for parallel processing
    df_samples <- tibble(xbar_a = as.vector(init_mat));
    for (n in seq_along(param_col_list)) { df_samples[[param_col_list[n]]] <- 
      as.vector(df_mat[[param_col_list[n]]]) }
    df_samples$diff_xbar_dm <- diff(xbars_dm)[1]
    
    # Compute desired statistics for each normalized samples
    # Control group 
    x_ctrl <- rnorm(n=df_samples$n_a[1], mean = 0, sd = 1)
    x_ctrl <- (x_ctrl - mean(x_ctrl))/sd(x_ctrl)
    x_ctrls <- t(matrix(rep(x_ctrl, dim(df_samples)[1]), ncol = dim(df_samples)[1]))
    x_ctrls <- sweep(sweep(x_ctrls, MARGIN=1, df_samples$sd_a, `*`), MARGIN=1, df_samples$xbar_a,'+')
    # Experiment group
    x_exp  <- rnorm(df_samples$n_b[2], mean = 0, sd = 1)
    x_exp <- (x_exp - mean(x_exp))/sd(x_exp)
    x_exps <- t(matrix(rep(x_exp, dim(df_samples)[1]), ncol = dim(df_samples)[1]))
    x_exps <- sweep(sweep(x_exps, MARGIN=1, df_samples$sd_b, `*`), MARGIN=1, df_samples$xbar_b,'+')
    
    # browser();
    
    # Calculate sample statistics
    if (stat_name == "mdm") { 
      df_samples[[stat_name]] <- row_mdm(x_ctrls, x_exps, conf.level = 1-alphas[1])
    }else if (stat_name == "rmdm") { 
      df_samples[[stat_name]] <- row_rmdm(x_ctrls, x_exps, conf.level = 1-alphas[1])
    } else {stop("process_cred_intervals: unsupported stat_name")}
    
    
    
    save(list = ls(all.names = TRUE), file = 
           "temp/process_cred_intervals.RData",envir = environment())
    # load(file = "temp/process_cred_intervals.RData")
    
    
    # Process parallel or serially
    if (is_parallel_proc) { 
      print("Starting parallel cluster...")
      n_nodes <- max(c(floor(detectCores()[1]*.9), 1))
      cl = makeCluster(n_nodes)
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
      
      ptm <- proc.time()
      df_lin2 <- 
        foreach(n = seq(1, n_xbars*n_sds, 1),
                .export = c("quant_cred_interval"),
                .combine = rbind,.packages = c("tidyverse")) %dopar% {
                  #calling a function
                  tempMatrix <- 
                    quant_cred_interval(df_samples[n,], stat_name = stat_name, n_mus = 50, n_sigmas = 52,
                                        p_thresh_mu = 1e-8, p_thresh_sigma = 1e-8, 
                                        method = method) 
                  
                  tempMatrix
                }
      est <- proc.time() - ptm
      
      # stopCluster(cl)
    } else {            # Process effect sizes serially
      # 
      
      df_list <- list()
      ptm <- proc.time()
      for (n in seq(1,n_xbars*n_sds,1)) {
        print(n)
        df_list[[n]] <- quant_cred_interval(df_samples[n,], stat_name = stat_name,
                                            n_mus = 25, n_sigmas = 25,
                                            p_thresh_mu = 1e-8, p_thresh_sigma = 1e-8, 
                                            method = method) 
      }
      sys_time <- proc.time() - ptm
      df_lin2 <- do.call("rbind", df_list)
    }
    
    
    
    
    save(list = ls(all.names = TRUE), file = "temp/process_cred_intervals2.RData",envir = environment())
    # load(file = "temp/process_cred_intervals2.RData")
    
    # Convert linear data frame to dataframe of matrices
    df_mat2 <- df_mat; col_list <- colnames(df_lin2)
    for (n in seq_along(col_list)) { 
      df_mat2[[col_list[n]]] <- matrix(df_lin2[[col_list[n]]], nrow = n_sds, 
                                       ncol = n_xbars, dimnames = dimnames)
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







# 
# 
# 
# 
# ## SCRATCH SPACE
# df_pops = tibble(mu_dm = seq(-5,5,.1), 
#                  sigma_dm = df_samples$sd_dm[1000],
#                  df = df_samples$n_a[1] + df_samples$n_b[1]-2 )
# df = quant_cred_interval(df_samples[1000,], df_pops, "mdm", "mu_dm") 
# 
# 
# 
# df_samples$
# # SCRATCH SPACE




