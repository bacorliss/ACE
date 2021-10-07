



#' QUantifies coverage error of the mdm (defined as how often abs(mu)>mdm_95 
#' from repeated samples)
#' Results computed in a grid with mu and sigma swept on each exis respectively

# Load required packages
#-------------------------------------------------------------------------------
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
p_load(binom)
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
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
# Script Parameters
fig_num = "4"
fig_path = file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE,recursive = TRUE)
n_samples <- 1e3
rand.seed <- 0
overwrite <- TRUE
is_parallel_proc <- TRUE







quant_cred_interval <-  function(df_sample, df_pops, sample_stat, pop_stat)  {
  #' @description Calculates the coverage error of x_expar, rx_expar, mdm, and rmdm
  #'
  #' @param df: a single row dataframe generated from agreement_contest
  #'  library, returned from generateExperiment_Data()
  #' @param n_samples: number of samples drawn to evaluate converage error
  #' @param n_obs: 
  #' @param enforce_grand_mean flag to force the grand mean (mean across all 
  #' samples and observations) to be exactly equal to the population mean. This 
  #' is used for debugging/ research purposes and not for actual simulations.
  #' @return null, exports figures to disk
  
  # save(list = ls(all.names = TRUE), file = "temp/quant_cred_interval.RData",envir = environment())
  # load(file = "temp/quant_cred_interval.RData")
  
  # Calculate probability that each pop param set can yield the specified sample mean
  # p_sample = dt(x = (df_sample$xbar_dm - df_pops$mu_dm)/df_pops$sd_dm, 
  #               df = df_sample$n_a + df_sample$n_b-2,
  #               ncp = df_pops$mu_dm/df_pops$sd_dm) * df_sample$diff_xbar_dm / df_sample$sd_dm
  
  
  # Integrate between tick marks in grid
  p_sample <-
    pt(q = (df_sample$xbar_dm - df_pops$mu_dm)/df_pops$sd_dm +
       0.5*df_sample$diff_xbar_dm / df_sample$sd_dm,
     df = df_sample$n_a + df_sample$n_b-2,
     ncp = df_pops$mu_dm/df_pops$sd_dm) -
    pt(q = (df_sample$xbar_dm - df_pops$mu_dm)/df_pops$sd_dm -
         0.5*df_sample$diff_xbar_dm / df_sample$sd_dm,
       df = df_sample$n_a + df_sample$n_b-2,
       ncp = df_pops$mu_dm/df_pops$sd_dm)
  
  # Calculate 
  p_pass  <- sum(p_sample[abs(df_pops$mu_dm) <= df_sample$stat], na.rm = TRUE)
  p_fail <- sum(p_sample[abs(df_pops$mu_dm) > df_sample$stat], na.rm = TRUE)
  
  df_out <- df_sample
  df_out[[paste("n_trials_", pop_stat, "_pass",sep="")]] = p_pass 
  df_out[[paste("n_trials_", pop_stat, "_fail",sep="")]] = p_fail
  df_out[[paste("credibility_rate_",pop_stat,"_lt_",sample_stat,sep="")]] = p_pass /(p_pass+p_fail)
  return(df_out)
}




quant_cred_intervals <- function(df_samples, df_pops, sample_stat, pop_stat)  {
  
  save(list = ls(all.names = TRUE), file = "temp/quant_cred_intervals.RData",envir = environment())
  # load(file = "temp/quant_cred_intervals.RData")
  
  df_list <- list()
  for (n in seq(1, dim(df_samples)[1],1)) {
    df_list[[n]] <- quant_cred_interval(df_samples[n,], df_pops, "mdm", "mu_dm") 
  }
  df_lin <- do.call("rbind", df_list)
  
  return(df_lin)
}





process_cred_intervals <-function(xbars_a, sds_a, n_a, xbars_b, sds_b, n_b, alphas = 0.05,
                             n_samples = n_samples, out_path,
                             overwrite = TRUE, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                             included_stats) {
  save(list = ls(all.names = TRUE), file = "temp/process_cred_intervals.RData",envir = environment())
  # load(file = "temp/process_cred_intervals.RData")
  dir.create(dirname(out_path),showWarnings = FALSE, recursive = TRUE)
  
  # Only compute requested statistics
  df_include = data.frame(matrix(rep(TRUE,length(included_stats)), nrow = 1, 
                                 dimnames = list(NULL, included_stats)))
  
  # For any pop param not equal in length to n_sims, expand
  n_xbars = max(length(xbars_b), length(xbars_a))
  n_sds = max(length(sds_a), length(sds_b))
  if (length(xbars_b)==1) {xbars_b = rep(xbars_b, n_xbars)}
  if (length(xbars_a)==1) {xbars_a = rep(xbars_a, n_xbars)} 
  if (length(sds_a)==1) {sds_a = rep(sds_a, n_sds)}
  if (length(sds_b)==1) {sds_b = rep(sds_b, n_sds)}
  
  # Caclulate difference in mean statistics
  xbar_dm = xbars_b - xbars_a
  sds_dm = sqrt(sds_a^2/n_a + sds_b^2/n_b) 
  
  # Store results to disk since calculations are significant
  set.seed(rand.seed)
  if (!file.exists(out_path) || overwrite) {
    
    # Matrix diff and error of mdm
    dimnames = list(sds_dm,xbar_dm)
    
    # Initialize matrix in data form
    param_col_list <- c("xbar_a", "sd_a","n_a", "xbar_b","sd_b", "n_b", 
                        "xbar_dm", "sd_dm", "rxbar_dm", "n_samples" ,"alpha")
    
    init_mat <- matrix(NA, nrow = n_sds, ncol = n_xbars, dimnames = dimnames)
    df_mat <- tibble(xbar_a = init_mat);
    for (n in seq_along(param_col_list)) {df_mat[[param_col_list[n]]] <- init_mat }
    
    for (r in seq(1, n_sds, by = 1)) {
      for (c in seq(1, n_xbars, by = 1)) {
        df_mat$xbar_a[r,c]     <- xbars_b[c]
        df_mat$sd_a[r,c]  <- sds_a[r]
        df_mat$n_a[r,c]  <- n_a
        
        df_mat$xbar_b[r,c]     <- xbars_a[c]
        df_mat$sd_b[r,c]  <- sds_b[r]
        df_mat$n_b[r,c]  <- n_b
        
        # Difference in means group
        df_mat$xbar_dm[r,c]     <- xbars_dm[c]
        df_mat$sd_dm[r,c]  <- sds_dm[r]
        df_mat$alpha[r,c]  <- alphas
        
        df_mat$rxbar_dm[r,c] <- df_mat$xbar_dm[r,c] / df_mat$xbar_a[r,c]
        df_mat$n_samples[r,c] = n_samples
        df_mat$n_a[r,c] = n_a
        df_mat$n_b[r,c] = n_b
        
      }
    } 
    
    # Linearize matrix dataframe for parallel processing
    df_samples <- tibble(xbar_a = as.vector(init_mat));
    for (n in seq_along(param_col_list)) { df_samples[[param_col_list[n]]] <- 
      as.vector(df_mat[[param_col_list[n]]]) }
    df_samples$diff_xbar_dm <- diff(xbar_dm)[1]
    
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
    
    
    # Calculate sample statistics
    if (!is.null(df_include$mdm)) { 
      df_samples$stat <- row_mdm(x_ctrls, x_exps, conf.level = 0.95)
    }
    if (!is.null(df_include$rmdm)) { 
      df_samples$stat <- row_rmdm(x_ctrls, x_exps, conf.level = 0.95)
    }
    

    # Define pop parameter space (must be much larger than sample parameter space)
    mu_dm_set <- seq(-6*max(as.vector(df_samples$stat)),
                     6*max(as.vector(df_samples$stat)), diff(xbar_dm)[1])
    sigma_dm_set <- seq(6*sds_dm[1], 6*sds_dm[length(sds_dm)], diff(sds_dm)[1])
    # Produce a grid for all mu_dm_set and sigma_dm_set combinations
    df_pops <- tibble(mu_dm = rep(mu_dm_set, length(sigma_dm_set)),
                     sd_dm = sigma_dm_set[rep(1:length(mu_dm_set), each=length(sigma_dm_set))],
                     df = df_samples$n_a[1] + df_samples$n_b[1]-2)
    n_mus = length(mu_dm_set)
    n_sigmas = length(sigma_dm_set)
    
    
    # Process parallel or serially
    if (is_parallel_proc) { 
      print("Starting parallel cluster...")
      n_nodes <- max(c(floor(detectCores()[1]*.9), 1))
      cl = makeCluster(n_nodes)
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
      
      save(list = ls(all.names = TRUE), file = 
             "temp/quant_cred_intervals.RData",envir = environment())
      
      ptm <- proc.time()
      df_lin2 <- 
        foreach(n = seq(1, n_xbars*n_sds, 1),
                .export = c("quant_cred_intervals", "quant_cred_interval"),
                .combine = rbind) %dopar% {
                  #calling a function
                  tempMatrix <- 
                    quant_cred_interval(df_samples[n,], df_pops, "mdm", "mu_dm") 
                  tempMatrix
                }
      est <- proc.time() - ptm
     
      # stopCluster(cl)
    } else {            # Process effect sizes serially
      # 
      
      df_list <- list()
      ptm <- proc.time()
      for (n in seq(1,n_sigmas*n_mus,1)) {
        print(n)
        df_list[[n]] <- quant_cred_interval(df_samples[n,], df_pops, "mdm", "mu_dm") 
      }
      sys_time <- proc.time() - ptm
      df_lin2 <- do.call("rbind", df_list)
    }
    
    
    save(list = ls(all.names = TRUE), file = "temp/quant_cred_intervals2.RData",envir = environment())
    # load(file = "temp/quant_cred_intervals2.RData")
    # browser();
    
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




# Row 1: 2D visualization of mdm difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#

n_obs = 6
xbars_dm <- seq(-2.5, 2.5, by = .1)
sds_dm <- seq(.01, 1, by = .05)

# Spread sigma_dm across sigma_a and sigma_b equally
sds_a = sds_dm/sqrt(2/n_obs)
sds_b = sds_a
n_samples <- 1e3


# Run simulations calculating error of mdm with mu and sigma swept
df_results <- 
  process_cred_intervals(xbars_a = 100, sds_a = sds_a, n_a = n_obs, 
                        xbars_b = 100 + xbars_dm, sds_b = sds_b, n_b = n_obs, alphas = 0.05,
                        n_samples = n_samples, out_path = paste(fig_path, "/mdm_cred_xbar_vs_s.rds",sep=""),
                        overwrite=overwrite, is_parallel_proc = TRUE, raw_error = TRUE, rel_error = FALSE,
                        included_stats = c("mdm"))








# 1A: Error rate of MDM < mu
#------------------------------------------------------------------------------#
# Convert from matrix to dataframe
df <- cbind(sigma = sds_dm, as_tibble(df_results$credibility_rate_mu_dm_lt_mdm)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
# grid_slopes <- slopes_by_rowcol(df_results$mean_mdm_error_rate, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu[DM])) + ylab(expression(sigma[DM])) +
  theme_classic(base_size=8) +
  scale_fill_gradientn(colors=c(rep("blue",18),"white", "red"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(0,1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste(fig_path, "/", fig_num, "_1a mdm error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 
