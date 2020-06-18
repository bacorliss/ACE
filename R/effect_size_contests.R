### Effect Size Contest



# TODO: for specifying offset distribution in generateData(), should the standard
# deviation for the control group be added by stds or variances?

library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
# library(gridExtra)
# library(grid)
library(tidyr)
library(cowplot)
library(dplyr)
# library(effsize)
library(boot)
# User defined libraries
source("R/row_effect_sizes.R")
source("R/mmd.R")

# Dictionary to keep track of variables tracking effect size metrics
effect_size_dict <- vector(mode="list", length=4)
names(effect_size_dict) <- c("prefix", "base", "suffix","label")
effect_size_dict[[1]] <- c("fract", "mean_diff")
effect_size_dict[[2]] <- c("xdbar", "rxdbar", "sd", "rsd", "bf", "p_value",
                           "cohen_d", "hedge_g", "glass_delta", "mmd",
                           "rmmd","nrand")
effect_size_dict[[3]] <- c("d2gtd1","2m1")
effect_size_dict[[4]] <- c("bar(x)", "r*bar(x)", "s", "r*s","Bf", "p",
                           "Cd", "Hg", "G*Delta", "delta[M]",
                           "r*delta[M]","Rand")


# default distribution for population parameters for Exp 1 {a,b}, Exp 2 {a,b}
generateExperiment_Data <- function(n_samples, n_obs, n_sims, rand.seed,
                                    # Control group pop. parameters
                                    mus_1a, sigmas_1a, 
                                    mus_2a, sigmas_2a,
                                    # Experiment group pop. parameters
                                    mus_1b = NA, sigmas_1b = NA, 
                                    mus_2b = NA, sigmas_2b = NA,
                                    # Difference distribution pop. parameters
                                    mus_1d, sigmas_1d, 
                                    mus_2d, sigmas_2d,
                                    label_dict = effect_size_dict) {
  #' Generate simulated experiment data for two experiments 
  #' 
  #' @description Generate simulated experiment data for two experiments with 
  #' distributions for mu and sigma specified as functions
  #' 
  #' @param n_samples number of samples (collection of observations), number of 
  #' times that a simulated experiment is repeated
  #' @param n_obs number of measurements in a sample/experiment
  #' @param n_sims number of simulations run, which are sets of experiments with
  #'  a single set of parameters for each (mu, sigma, n_obs)
  #' @param rand.seed seed number of generating consistent random numbers
  #' @param mus1 vector specifying distribution for population mean for
  #'  experiment 1
  #' @param sigmas1 vector specifying distribution for population standard deviation for
  #'  experiment 1
  #' @param mus2 vector specifying distribution for population mean for
  #'  experiment 2
  #' @param sigmas2 vector specifying distribution for population mean for
  #'  experiment 12
  
  # Basic input parameters for each simulation round
  set.seed(rand.seed +1)
  df = tibble(
    # Sampling data
    n_obs = n_obs, n_samples = n_samples,
    # Control group a for exp 1 and 2
    mu_1a = mus_1a, mu_2a = mus_2a,
    sigma_1a = sigmas_1a, sigma_2a = sigmas_2a,
  )
  
 
  # Fill in experiment group and/or calculate difference distribution
  if (!any(is.na(mus_1b))) { 
    ##  Experiment group valid
    # Fill in experiment group
    # Experiment group b for exp 1 and 2
    df$mu_1b <- mus_1b
    df$sigma_1b <- sigmas_1b
    df$mu_2b <- mus_2b
    df$sigma_2b <- sigmas_2b

    # Difference distribution    
    df$mu_1d <- df$mu_1b - df$mu_1a
    df$mu_2d <- df$mu_2b - df$mu_2a
    df$sigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2)
    df$sigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) 


    
  } else {
    ##  Experiment group invalid
    
    # Calculate from difference distribution
    
    # Experiment group mean based on control and offset
    df$mu_1b <- mus_1a + mus_1d
    df$mu_2b <- mus_2a + mus_2d
    # Variance of difference in experimental
    df$sigma_1b <- sqrt(df$sigma_1a^2 + sigmas_1d^2)
    df$sigma_2b <- sqrt(df$sigma_2a^2 + sigmas_2d^2)
    
    # Calculate difference distribution parameters (not mean difference)
    df$mu_1d <- mus_1d
    df$mu_2d <- mus_2d
    # Variance of difference in means based on control and offset
    df$sigma_1d <- sigmas_1d
    df$sigma_2d <- sigmas_2d
    
  
    
    #browser();
  }
  
  # Calculate parameter of difference in means distribution
  df$mu_1md <- mus_1d
  df$mu_2md <- mus_2d
  df$sigma_1md <- df$sigma_1d /sqrt(n_obs)
  df$sigma_2md <- df$sigma_2d /sqrt(n_obs)
  
  # Calculate ratio of sigma_md/mu_md to determine how close D is close to zero,
  # and how absolute value folding will effect distribution.
  df$mu_ov_sigma_1md <- df$mu_1md / df$sigma_1md
  df$mu_ov_sigma_2md <- df$mu_2md / df$sigma_2md
  
  
  # Is: Exp2 mu[d] > Exp1 mu[d]
  df$is_mud_md2gtmd1 <-  abs(df$mu_2md) > abs(df$mu_1md)
  # Statistics of difference of means distribution 
  df$rmu_1md <- df$mu_1md/df$mu_1a
  df$rmu_2md <- df$mu_2md/df$mu_2a
  # Is: Exp2 relative_mu[d] > Exp1 relative_mu[d]
  df$is_rmud_md2gtmd1 <-  abs(df$rmu_2md) > abs(df$rmu_1md)
  
  # Is: pop_std2 > pop_std1 (TRUE)
  df$is_sigma_md2gtmd1 <-  df$sigma_2md > df$sigma_1md
  # Relative change variance of the mean compared to difference in means
  #   Realtive variance across scales
  # Is: rel_pop_std2 > rel_pop_std1, a/k/a CV1 > CV2
  df$rsigma_1md <- df$sigma_1md / abs(df$mu_1a + df$mu_1d/2)
  df$rsigma_2md <- df$sigma_2md / abs(df$mu_2a + df$mu_2d/2)
  df$is_rsigma_md2gtmd1 <-  df$rsigma_2md > df$rsigma_1md

  
  # Diff:  pop_mean2 - pop_mean1
  df$mean_mud_d2md1 <- df$mu_2md - df$mu_1md
  df$mean_rmud_d2md1 <- (df$mu_2md/df$mu_2a) - (df$mu_1md/df$mu_1a)
  # Diff:  pop_std2 - pop_std1
  df$mean_sigma_d2md1 <- df$sigma_2md - df$sigma_1md
  df$mean_rsigma_d2md1 <- df$sigma_2md/df$mu_2a - df$sigma_1md/df$mu_1a
  
  # Append columns for effect sizes, since multiple columns are used to analyze
  # each effect size, a dictionary of prefix, base, and suffix variable names 
  # are used.
  df[ paste(effect_size_dict$prefix[1], effect_size_dict$base,
            effect_size_dict$suffix[1], sep="_") ] <- rep(NaN,n_sims)
  df[ paste(effect_size_dict$prefix[2], effect_size_dict$base,
            effect_size_dict$suffix[2], sep="_") ] <- rep(NaN,n_sims)
  
  # browser()
  return(df)
}

quantify_esize_simulations <- function(df, overwrite = FALSE,
                                out_path = "temp/", data_file_name,
                                rand.seed = 0, include_bf = FALSE) {
  #' Simulate experiments generated from generateExperiment_Data() and calculates
  #'  various effect sizes
  #' 
  #' @description Given a data frame of parameters and input data for each 
  #' experiment group, quantifies data with various effect size metrics.
  #' 
  #' @param df input data frame generated from generateExperiment_Data() that 
  #' holds parameters for simulated data, along with initialized fields to store
  #'  effect size metrics
  #' @param x_a base input data from group a (gets transformed into actual input
  #'  data with parameters within df)
  #' @param x_b base input data from group b (gets transformed into actual input
  #'  data with parameters within df)
  #' @param out_path file path to save results to disk
  #' @param overwrite if results file already exists in out_path, skip 
  #' calculation and load from disk
  n_sims = dim(df)[1]
  
  # Only perform simulations if results not saved to disk
  if (!file.exists(paste(out_path,data_file_name,sep="")) | overwrite) {
    for (n in seq(1,n_sims,1)) {
      set.seed(rand.seed+n)
      # Transform simulated samples with normal parameteres (mean and std) for 
      # current round of simulation
      
      # Use Exp 1 and 2 coefficients to generate data from normalized base data
      x_1a = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_1a[n], 
                          sd = df$sigma_1a[n]), nrow = df$n_samples[n], 
                    ncol = df$n_obs[n])
      x_1b = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_1b[n], 
                          sd = df$sigma_1b[n]), nrow = df$n_samples[n], 
                    ncol = df$n_obs[n])
      
      x_2a = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_2a[n], 
                          sd = df$sigma_2a[n]), nrow = df$n_samples[n], 
                    ncol = df$n_obs[n])
      x_2b = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_2b[n], 
                          sd = df$sigma_2b[n]), nrow = df$n_samples[n], 
                    ncol = df$n_obs[n])
      
      # Sample estimate of difference in means
      xbar_1d = rowMeans(x_1b) - rowMeans(x_1a)
      xbar_2d = rowMeans(x_2b) - rowMeans(x_2a)
      
      # Sample estimate of standard deviation of difference in means
      s_1md = sqrt(rowSds(x_1a)^2/df$n_obs[n] + rowSds(x_1b)^2/df$n_obs[n])
      s_2md = sqrt(rowSds(x_2a)^2/df$n_obs[n] + rowSds(x_2b)^2/df$n_obs[n])
      
      # Basic Summary statistical comparisons
      # Means
      df$fract_xdbar_d2gtd1[n] = sum(abs(xbar_2d) > abs(xbar_1d))/df$n_samples[n]
      df$mean_diff_xdbar_2m1[n]  = mean(abs(xbar_2d) - abs(xbar_1d))
      
      # Stds
      df$fract_sd_d2gtd1[n]  = sum(  s_2md > s_1md) / df$n_samples[n]
      df$mean_diff_sd_2m1[n] = mean( s_2md - s_1md)
      
      # Rel Means: mean divided by control mean
      diff_rxdbar = abs(xbar_2d/rowMeans(x_2a)) - abs(xbar_1d/rowMeans(x_1a))
      df$fract_rxdbar_d2gtd1[n] = sum( diff_rxdbar > 0) / df$n_samples[n]
      df$mean_diff_rxdbar_2m1[n]  =  mean(diff_rxdbar)
      
      # Rel STDs: sd divided by control mean
      diff_rsd = s_2md / (rowMeans(x_2a) + xbar_2d/2) -
        s_1md / (rowMeans(x_1a) + xbar_1d/2)
      df$fract_rsd_d2gtd1[n] = sum( diff_rsd > 0) / df$n_samples[n]
      df$mean_diff_rsd_2m1[n]  = mean(diff_rsd)
      
      if (df$fract_rsd_d2gtd1[n] > 1) {browser();}
      
      if (include_bf) {
        # Bayes factor
        bf_1 = row_ttestBF(x_1a, x_1b, parallel = TRUE, paired = FALSE)
        bf_2 = row_ttestBF(x_2a, x_2b, parallel = TRUE, paired = FALSE)
        diff_abs_bf <- abs(bf_2) - abs(bf_1)
        df$fract_bf_d2gtd1[n] = sum(diff_abs_bf > 0) /
          df$n_samples[n]
        df$mean_diff_bf_2m1[n] = mean(diff_abs_bf)
      } else {
        df$fract_bf_d2gtd1[n] = 0
        df$mean_diff_bf_2m1[n] = 0
      }
      
      
      # Pvalue
      # The more equal experiment will have a larger p-value
      z_score_1 <- rowzScore(x_1b, x_1a)
      z_score_2 <- rowzScore(x_2b, x_2a)
      p_value_1 = 2*pnorm(-abs(z_score_1))
      p_value_2 = 2*pnorm(-abs(z_score_2))
      diff_p_value <-  p_value_2 - p_value_1
      # Must switch signs
      df$fract_p_value_d2gtd1[n] = sum(-diff_p_value > 0) / df$n_samples[n]
      df$mean_diff_p_value_2m1 = mean(diff_p_value)
      
      # Delta Family of effect size
      # Cohens D
      diff_cohen_d = abs(rowCohenD(x_2a, x_2b)) - abs(rowCohenD(x_1a, x_1b))
      df$fract_cohen_d_d2gtd1[n] = sum(diff_cohen_d > 0) / df$n_samples[n]
      df$mean_diff_cohen_d_2m1[n] =  mean(diff_cohen_d)
      # Glass delta
      diff_glass_delta = abs(rowGlassDelta(x_2a, x_2b)) - abs(rowGlassDelta(x_1a, x_1b))
      df$fract_glass_delta_d2gtd1[n] = sum(diff_glass_delta > 0) / df$n_samples[n]
      df$mean_diff_glass_delta_2m1[n] =  mean(diff_glass_delta)
      # Hedges G
      diff_hedge_g = abs(rowHedgeG(x_2a, x_2b)) - abs(rowHedgeG(x_1a, x_1b))
      df$fract_hedge_g_d2gtd1[n] = sum(diff_hedge_g > 0) / df$n_samples[n]
      df$mean_diff_hedge_g_2m1[n] =  mean(diff_hedge_g)
      
      # Most Mean Diff
      mmd_1 = row_mmd(x_1a, x_1b, paired = FALSE)
      mmd_2 = row_mmd(x_2a, x_2b, paired = FALSE)
      diff_most_mean_diff = mmd_2 - mmd_1
      df$fract_mmd_d2gtd1[n] = sum(diff_most_mean_diff > 0) / df$n_samples[n]
      df$mean_diff_mmd_2m1[n] = mean(diff_most_mean_diff)
      
      # Relative Most Mean Diff
      diff_rmmd = mmd_2 / rowMeans(x_2a) -
        mmd_1 / rowMeans(x_1a)
      df$fract_rmmd_d2gtd1[n] = sum(diff_rmmd > 0) / df$n_samples[n]
      df$mean_diff_rmmd_2m1[n] = mean(diff_rmmd)
      
      
      # Random
      diff_nrand = rowMeans(matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = 0, sd = 1), 
                                   nrow = df$n_samples[n], ncol = df$n_obs[n]))
      df$fract_nrand_d2gtd1[n] = sum(diff_nrand > 0 ) / df$n_samples[n]
      df$mean_diff_nrand_2m1[n] = mean(diff_nrand)

      # Pearson Correlation
      # R Squared
      # Biserial Correlation https://rpubs.com/juanhklopper/biserial_correlation
    }
    # Save dataframed results to a file
    saveRDS(df, file = paste(out_path,data_file_name,sep=""))
  } else {
    # Restore the dataframed results from disk
    df <- readRDS(file = paste(out_path,data_file_name,sep=""))
    
  }
  # browser()
  return(df)
}



tidy_esize_simulations <- function (df, gt_colname, var_suffix, long_format = TRUE,
                                   ref_colname = NULL) {
  #' Normalize a subset of variables in df and then flatten data frame
  #' 
  #' @description Normalize a subset of variables in df, subract/compare to groundtruth, and 
  #' then subtract out the reference value if included
  #' 
  #' @param df
  #' @param gt_colname
  #' @param var_suffix
  #' @param long_format
  #' @param ref_colname
  
  # Check if gt_colname exists
  if (length(grep(gt_colname, names(df))) == 0) stop("gt_colname does not exist in dataframe")
  
  # Get list of all variables that match the variable suffix string
  matched_vars <- grep(var_suffix,colnames(df), perl=TRUE, value=TRUE)
  
  # Initialize new zeroed df with only selected vars
  data = matrix(0, dim(df)[1], length(matched_vars))
  colnames(data) <- matched_vars
  df_gt_sub <- as_tibble(data)
  
  # Fill in resulting values depending on variable type
  # Frac: fraction of metrics of samples that satisfy some condition: to process,
  #        subtract these values from ground truth.
  # mean_diff: mean value of difference between groups:  to process,
  #        test for equaivalent logical values
  for (n in seq(1,length(matched_vars), by = 1)) {
    if (var_suffix == "fract") {
      # If var_suffix is "fract.*", for cases where GT vector is TRUE, compute (1-x)
      df_gt_sub[[matched_vars[n]]][ df[[gt_colname]]] <- 
        1 - df[[matched_vars[n]]][ df[[gt_colname]]]
      df_gt_sub[[matched_vars[n]]][!df[[gt_colname]]] <-   
        df[[matched_vars[n]]][!df[[gt_colname]]]
    } else {
      # if var_suffix is "mean_diff.*", subtract from GT
      df_gt_sub[matched_vars[n]] = df[matched_vars[n]] - gt_vector
    }
  }
  
  # If reference value included, subtract in paired fashion (reference value 
  # must be included in matched_variable list, with groundtruth subtract from it)
  df_ref_sub <- df_gt_sub
  if (!is.null(ref_colname)) {
    for (n in seq(1,length(matched_vars), by = 1)) {
      df_ref_sub[matched_vars[n]] <- df_gt_sub[matched_vars[n]] - 
        df_gt_sub[ref_colname]
    }
  }
  
  # Flatten df_ref_sub into df_tidy (wide format to long format data frame)
  if (long_format) {
    df_tidy <- df_ref_sub %>% gather(matched_vars,factor_key=TRUE)
    names(df_tidy)[names(df_tidy) == "matched_vars"] <- "name"
  } else df_tidy <- df_ref_sub
  
  # browser()
  return(df_tidy)
}

pretty_esize_levels<- function(df,base_names, pretty_names, var_suffix) {
  #' Convert long levels names to shorter ones for pretty print in plots
  #' 
  #' @description 
  #' @param 
  #' @param 
  #' @param 
  #' 
  
  # Get current names for levels and the basenames for them to be matched against
  orig_levels <- levels(df$name)
  # base_names <- paste(var_suffix,"_", base_names, "_", sep="") 
  
  # Find index for each basename
  ind <- rep(0, length(base_names))
  for (n in seq(1,length(base_names))) {
    new_levels_ind =1
    ind[n] <- which(regexpr(base_names[n],orig_levels)>0)
  }
  # Return new level names
  new_levels <- effect_size_dict$label[ind]
  
  # df_new_levels <- df
  # levels(df_new_levels$name) <- new_levels
  levels(df$name) <- new_levels
  
  return(df)
  
}


plot_population_params <- function(df_init, gt_colname,fig_name,y_ax_str){
  
  # Plot reference ground truth success rate (Exp 1 < Exp 2)
  # Check to see overall ground truth true rate
  binom_p <- prop.test(sum(df_init[[gt_colname]]), dim(df_init)[1], conf.level=0.95, correct = FALSE)
  p <- ggplot(tibble(x=as.factor(1),y=binom_p$estimate), aes(x=x,  y=y)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey") +
    geom_linerange(aes(ymin = binom_p$conf.int[1], ymax = binom_p$conf.int[2]), size = 0.5) +
    geom_point(size=1,fill="white", shape=1) + 
    ylab("Fract Exp 1 < Exp 2") +
    xlab( parse(text = y_ax_str)) +
    theme_classic(base_size = 8) +
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  p
  save_plot(paste("figure/", 'gt_',fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = .75, dpi = 600)
  
  
  # Plot histogram of mu[D]/sigma[D] to demonstrate how far from zero D is  
  df <-tibble(group = as.factor(c(rep(1,dim(df_init)[1]),rep(2,dim(df_init)[1]))),
              mu_ov_sigma = abs(c(df_init$mu_1md/df_init$sigma_1md,
                                  df_init$mu_2md/df_init$sigma_2md)))
  p <- ggplot(df, aes(x = mu_ov_sigma, y = mu_ov_sigma, fill = group)) +
    geom_histogram(aes(y=stat(count / sum(count))), position="identity", 
                   alpha=0.25, bins = 15) +
    xlab( expression(abs(~mu[D]*phantom(.))*phantom(.)/phantom(.)*sigma[D])) +
    ylab( "Freq.") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  #print(p)
  save_plot(paste("figure/", 'mu_ov_sigma_',fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 1.5, dpi = 600)
  
  
}

plot_esize_simulations <- function(df_pretty, fig_name, y_ax_str) {
  # Calculate confidence interval
  df_result <- df_pretty %>%   
    group_by(name) %>% 
    summarize(mean = mean(value), bs_ci_mean_str = toString(
      boot.ci(boot(value, function(x, ind)  mean(x[ind]), R = 10000),
              conf = 1-(0.05/length(levels(df_pretty$name))), type = "basic" )$
        basic[c(4,5)]))
  df_result$bs_ci_mean_lower <- sapply(strsplit(df_result$bs_ci_mean_str,","), 
                                       function(x) as.numeric(x[1]))
  df_result$bs_ci_mean_upper <- sapply(strsplit(df_result$bs_ci_mean_str,","), 
                                       function(x) as.numeric(x[2]))
  # If groups have no variance in mean, then boot returns NaNs, replace with mean
  df_result$bs_ci_mean_lower[is.na(df_result$bs_ci_mean_lower)] <- 
    df_result$mean[is.na(df_result$bs_ci_mean_lower)]
  df_result$bs_ci_mean_upper[is.na(df_result$bs_ci_mean_upper)] <- 
    df_result$mean[is.na(df_result$bs_ci_mean_upper)]
  
  df_result$is_mean_0.5 <- 0.5 >= df_result$bs_ci_mean_lower &
    0.5 <= df_result$bs_ci_mean_upper
  ci_range <-  c(min(df_result$bs_ci_mean_lower), max(df_result$bs_ci_mean_upper))
  
  # Labels of statistical significance for each group
  sig_labels = rep("",length(levels(df_pretty$name)))
  sig_colors = rep("black",length(levels(df_pretty$name)))
  sig_sizes = rep(4,length(levels(df_pretty$name)))
  siff_vjust = rep(0,length(levels(df_pretty$name)))
  # Set less than random to blue and -
  sig_labels[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  "-"
  sig_colors[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  rgb(47, 74, 71,maxColorValue = 255)
  sig_sizes[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  5
  siff_vjust[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  .017
  # Set greater than random to red and +
  sig_labels[df_result$bs_ci_mean_lower>0.5 & df_result$bs_ci_mean_upper>0.5] =  "+"
  sig_colors[df_result$bs_ci_mean_lower>0.5 & df_result$bs_ci_mean_upper>0.5] =  rgb(255, 0, 0,maxColorValue = 255)
  
  
  # Basic violin plot
  p <- ggplot(df_result, aes(x=name,  y=mean, group=name)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey") +
    geom_linerange(aes(ymin = bs_ci_mean_lower, ymax = bs_ci_mean_upper), size = 0.5) +
    geom_point(size=1,fill="white", shape = 1) + 
    xlab("Effect Size Metric") +
    ylab(parse(text=paste("Error~Rate~Lower~phantom(.)*", y_ax_str))) +
    scale_x_discrete(labels = parse(text = levels(df_pretty$name))) +
    expand_limits(y = extend_max_lim(ci_range, 0.2)) +
    geom_text(y = extend_max_lim(ci_range, 0.2)+siff_vjust, aes(label = sig_labels), 
              color = sig_colors, size = sig_sizes, vjust=0.5, hjust=0.5) +
    theme_classic() +  theme(text = element_text(size = 8))
  print(p)
  save_plot(paste("figure/", fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 3.25, dpi = 600)
   #browser()

  return(df_result)
  
}

process_esize_simulations <- function(df_init, gt_colname, y_ax_str, out_path="temp/",
                                      fig_name,var_suffix = "fract") {
  
  # Display ground truth fraction of E2>E1
  print(sprintf("%s (TRUE): %i", gt_colname, sum(df_init[[gt_colname]])))
  
  # Plot data about population params
  plot_population_params(df_init, gt_colname,fig_name,y_ax_str)
    
    
  # Quantify effect sizes in untidy matrix
  df_es <- quantify_esize_simulations(df = df_init,overwrite = TRUE, out_path = out_path,
                                      data_file_name = paste(fig_name,".rds",sep = ""))
  
  # Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
  df_tidy <- tidy_esize_simulations(df = df_es, gt_colname = gt_colname,
                                   var_suffix = var_suffix,long_format = TRUE,
                                   ref_colname = NULL)
  df_pretty <- 
    pretty_esize_levels(df = df_tidy, base_names = paste(var_suffix,"_",
                                                         effect_size_dict$base, "_", sep=""),
                                         pretty_names = effect_size_dict$label, 
                                         var_suffix = var_suffix)
  # browser()
  # Plot effect size results
  df_plotted <- plot_esize_simulations(df = df_pretty, fig_name = fig_name, y_ax_str = y_ax_str)
  
  
  all_dfs <- vector(mode="list", length=4)
  names(all_dfs) <- c("df_es", "df_tidy", "df_pretty", "df_plotted")
  all_dfs[[1]] <- df_es; all_dfs[[2]] <- df_tidy; 
  all_dfs[[3]] <- df_pretty; all_dfs[[4]] <- df_plotted; 
  return(all_dfs)
}