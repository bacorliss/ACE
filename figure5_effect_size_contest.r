

### Figure 2: MMD outperforms previous statistics of effect size.
#
# Metrics of sample effect size
#
##############################################
#  Unstandardized measures of effect size
#   * difference of sample means
#   * relative difference of means
#   * difference of sample variances
#  Standardized measures of effect size
#   * cohen's d
#   * hedges g
#   * glasses delta
#   * MHD


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

# User defined libraries
source("R/mmd.R")

## User defined functions
# Calculate variance by row like rowMeans or rowSums
rowVars <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)}
rowSds  <- function(x, ...) sqrt(rowVars(x))
# Pooled standard deviation, assume to matricies, 1 sample per row
rowSdPooled <- function(m1,m2) sqrt( (rowSds(m1)^2 + rowSds(m2)^2 )/2 )
# Weighted pooled standard deviation, assume to matricies, 1 sample per row
rowSdPooledWeighted <- 
  function(m1,m2) sqrt( ( (dim(m1)[2] - 1) * rowVars(m1)   + 
                            (dim(m2)[2] - 1) * rowVars(m2) ) /
                          (dim(m1)[2] + dim(m2)[2] - 2    ) )


# Effect Size Statistics Functions
#
#------------------------------------------------------------------------------
#
#   d = (M2 - M1)/s_pool
rowCohenD <- function (m1,m2) ( rowMeans(m1) - rowMeans(m2)) / rowSdPooled(m1,m2)
#   d = (M2 - M1)/s_pool
rowHedgeG <- function (m1,m2) ( rowMeans(m1) - rowMeans(m2)) / rowSdPooledWeighted(m1,m2)
#   d = (M2 - M1)/s_1
rowGlassDelta <- function (m1,m2) ( rowMeans(m2) - rowMeans(m1)) / rowSds(m1) 

rowTScore <- function(m1, m2) {
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  n  <- n1+n2 
  tstat<- sqrt(n1*n2/n) * ( rowMeans(m1) - rowMeans(m2)) /
    sqrt( (n1 - 1)/(n - 2)*rowVars(m1) + (n2 - 1)/(n - 2)*rowVars(m2) ) 
}

# Simulation Functions
#
#-------------------------------------------------------------------------------
#
# Code structure: 
# 1. two groups of untransformed sample sets are generated (x_a,x_b)
# 2. Transform coefficients for experiment 1 and 2 are generated
# 3. For each simulation, coefficients are applied
#     1. x_a [exp_1_coeffs]-> x1a
#     2. x_b [exp_1_coeffs]-> x1b
#     1. x_a [exp_2_coeffs]-> x2a
#     2. x_b [exp_2_coeffs]-> x2b
#     
#                      x_a                  x_b
#                   |       \           /        |
#                   |         \       /          |
#                   |           \   /            |
#                   |            / \             |
#                   |          /     \           |
#                ______________    ______________
#               | exp_1_coeffs |  | exp_2 coeffs |
#               |______________|  |______________| 
#                   |       /           \      | 
#                 x_1a    x_1b        x_2a    x_2b
#                    \    /               \   /
#                     \  /                 \ /
#                   metric_1            metric_2
#
#
# Results are store in a dataframe, where each row is a simulation with summary 
#     data from a large collection of samples. For each metric
#
# 1. fraction of samples where exp 2 has a greater metric value than exp 1
# 2. mean difference between metric values of [exp 2] - [exp 1]


# Dictionary to keep track of variables tracking effect size metrics
effect_size_dict <- vector(mode="list", length=4)
names(effect_size_dict) <- c("prefix", "base", "suffix","label")
effect_size_dict[[1]] <- c("fract", "mean_diff")
effect_size_dict[[2]] <- c("xdbar", "sd", "rxdbar","rsd", "t_score", "p_value",
                           "cohen_d", "hedge_g", "glass_delta", "mmd",
                           "rmmd","nrand")
effect_size_dict[[3]] <- c("d2gtd1","2m1")
effect_size_dict[[4]] <- c("bar(x)", "s", "r*bar(x)","r*s","t", "p",
                           "Cd", "Hg", "G*Delta", "delta[M]",
                           "r*delta[M]","N(0,1)")

# default distribution for population paramters for Exp 1 {a,b}, Exp 2 {a,b}



generateExperiment_Data <- function(n_samples, n_obs, n_sims, rand.seed,
                                    # Control group pop. parameters
                                    mus_1a  = rep(0,n_sims), 
                                    sigmas_1a = rep(1,n_sims), 
                                    mus_2a  = rep(0,n_sims), 
                                    sigmas_2a = rep(1,n_sims),
                                    # Experiment group pop. parameters
                                    mus_1b  = runif(n_sims, -0.5, 0.5), 
                                    sigmas_1b = runif(n_sims,  0.1, 1), 
                                    mus_2b  = runif(n_sims, -0.5, 0.5), 
                                    sigmas_2b = runif(n_sims,  0.1, 1),
                                    label_dict = effect_size_dict) {
  #' Generate simulated experiment data for two experiments 
  #' 
  #' @description Generate simualted experiment data for two experiments with 
  #' distributions for mu and sigma specified as functions
  #' 
  #' @param n_samples numer of samples (collection of observations), number of 
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
    # Control group a for exp 1 and 2
    mu_1a = mus_1a, sigma_1a = sigmas_1a,
    mu_2a = mus_2a, sigma_2a = sigmas_2a,
    # Experiment group b for exp 1 and 2
    mu_1b = mus_1b, sigma_1b = sigmas_1b,
    mu_2b = mus_2b, sigma_2b = sigmas_2b,
    # Sampling data
    n_obs = n_obs, n_samples = n_samples
  )
  # Computed statistics of difference of means distribution (mu2 - mu1)
  # Is: Exp2 mu[d] > Exp1 mu[d]
  df$is_mud_d2gtd1 <-  (df$mu_2b - df$mu_2a) > (df$mu_1b - df$mu_1a)
  # Is: Exp2 relative_mu[d] > Exp1 relative_mu[d]
  df$is_rmud_d2gtd1 <-  (df$mu_2b - df$mu_2a)/df$mu_2a > (df$mu_1b - df$mu_1a)/df$mu_1a
  # Diff:  pop_mean2 - pop_mean1
  df$mean_mud_d2md1 = (df$mu_2b - df$mu_2a) - (df$mu_1b - df$mu_1a)
  df$mean_rmud_d2md1 = (df$mu_2b - df$mu_2a)/df$mu_2a - (df$mu_1b - df$mu_1a)/df$mu_1a
  # Is: pop_std2 > pop_std1 (TRUE)
  df$is_sigmad_d2gtd1 =  (df$sigma_2b - df$sigma_2a) > (df$sigma_1b - df$sigma_1a)
  # Is: rel_pop_std2 > rel_pop_std1 (TRUE), AKA coefficient of variation
  df$is_rsigmad_d2gtd1 =  (df$sigma_2b - df$sigma_2a)/df$mu_2a > 
    (df$sigma_1b - df$sigma_1a)/df$mu_1a
  # Diff:  pop_std2 - pop_std1
  df$mean_sigmad_d2md1 = (df$sigma_2b - df$sigma_2a) - (df$sigma_1b - df$sigma_1a)
  df$mean_rsigmad_d2md1 = (df$sigma_2b - df$sigma_2a)/df$mu_2a -
    (df$sigma_1b - df$sigma_1a)/df$mu_1a
  
  # Append columns for effect sizes, since multiple columns are used to analyze
  # each effect size, a dictionary of prefix, base, and suffix variable names 
  # are used.
  df[ paste(effect_size_dict$prefix[1], effect_size_dict$base,
            effect_size_dict$suffix[1], sep="_") ] <- rep(NaN,n_sims)
  df[ paste(effect_size_dict$prefix[2], effect_size_dict$base,
            effect_size_dict$suffix[2], sep="_") ] <- rep(NaN,n_sims)
  return(df)
}

quantifyEffectSizes <- function(df, overwrite = FALSE,
                                out_path = "temp/quanitfyEffectSizes.rds", 
                                rand.seed = 0) {
  #' Simulate experiments generated from generateExperiment_Data() and calcualte
  #'  various effecct sizes
  #' 
  #' @description Given a data frame of parameters and input data for each 
  #' experiment group, quantifies data with various effect size metrics.
  #' 
  #' @param df input data frame generated from generateExperiment_Data() that 
  #' holds parameters for simualted data, along with initialized fields to store
  #'  effect size metricss
  #' @param x_a base input data from group a (gets transformed into actual input
  #'  data with parameters within df)
  #' @param x_b base input data from group b (gets transformed into actual input
  #'  data with parameters within df)
  #' @param out_path file path to save results to disk
  #' @param overwrite if results file already exists in out_path, skip 
  #' calculation and load from disk
  n_sims = dim(df)[1]
  
  # Only perform simulations if results not saved to disk
  if (!file.exists(out_path) | overwrite) {
    for (n in seq(1,n_sims,1)) {
      set.seed(rand.seed+n)
      # Transform simulated samples with normal parameteres (mean and std) for 
      # current round of simulation
      
      # Use Exp 1 and 2 coefficients to generate data from normalized base data
      x_1a = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_1a[n], sd = df$sigma_1a[n]),
                    nrow = df$n_samples[n], ncol = df$n_obs[n])
      x_1b = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_1b[n], sd = df$sigma_1b[n]),
                    nrow = df$n_samples[n], ncol = df$n_obs[n])
      
      x_2a = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_2a[n], sd = df$sigma_2a[n]),
                    nrow = df$n_samples[n], ncol = df$n_obs[n])
      x_2b = matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = df$mu_2b[n], sd = df$sigma_2b[n]),
                    df$n_samples[n], df$n_obs[n])
      
      # Calculate difference of basic statistics
      # Metric[B] - Metric[A]
      xdbar_1 = abs(rowMeans(x_1b) - rowMeans(x_1a))
      xdbar_2 = abs(rowMeans(x_2b) - rowMeans(x_2a))
      
      sdd_1 = rowSds(x_1b) - rowSds(x_1a)
      sdd_2 = rowSds(x_2b) - rowSds(x_2a)
      
      # Basic Summary statistical comparisons
      # Means
      df$fract_xdbar_d2gtd1[n] = sum(xdbar_2 > xdbar_1)/df$n_samples[n]
      df$mean_diff_xdbar_2m1[n]  = mean(xdbar_2 - xdbar_1)
      
      # Stds
      df$fract_sd_d2gtd1[n] = sum(sdd_2 > sdd_1)/df$n_samples[n]
      df$mean_diff_sd_2m1[n]  = mean(sdd_2 - sdd_1)
      
      # Rel Means: mean divided by control mean
      diff_rxdbar = xdbar_2/rowMeans(x_2a) - xdbar_1/rowMeans(x_1a)
      df$fract_rxdbar_d2gtd1[n] = sum( diff_rxdbar > 0) / df$n_samples[n]
      df$mean_diff_rxdbar_2m1[n]  =  mean(diff_rxdbar)
      
      
      # Rel STDs: sd divided by control mean
      diff_rsd = sdd_2/rowMeans(x_2a) - sdd_1/rowMeans(x_1a)
      df$fract_rsd_d2gtd1[n] = sum( diff_rsd > 0) / df$n_samples[n]
      df$mean_diff_rsd_2m1[n]  = mean(diff_rsd)
      
      
      # t score
      t_score_1 <- rowTScore(x_1b, x_1a)
      t_score_2 <- rowTScore(x_2b, x_2a)
      df$fract_t_score_d2gtd1[n] = sum( abs(t_score_2) - abs(t_score_1) > 0) / df$n_samples[n]
      df$mean_diff_t_score_2m1[n] = mean( (t_score_2 - t_score_1))
      
      # Pvalue
      diff_p_value <- pt(t_score_2, df = pmin(df$n_obs[n], df$n_obs[n]) - 1) - 
        pt(t_score_1, df = pmin(df$n_obs[n], df$n_obs[n]) - 1)
      df$fract_p_value_d2gtd1[n] = sum(diff_p_value > 0) / df$n_samples[n]
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
      most_mean_diff_1 = sapply(1:df$n_samples[n], function(i) 
        mmd_normal(x_1a[i,], x_1b[i,],paired = FALSE))
      most_mean_diff_2 = sapply(1:df$n_samples[n], function(i) 
        mmd_normal(x_2a[i,], x_2b[i,], paired = FALSE))
      diff_most_mean_diff = most_mean_diff_2 - most_mean_diff_1
      df$fract_mmd_d2gtd1[n] = sum(diff_most_mean_diff > 0) / df$n_samples[n]
      df$mean_diff_mmd_2m1[n] = mean(diff_most_mean_diff)
      
      # Relative Most Mean Diff
      diff_rmmd = most_mean_diff_2/rowMeans(x_2a) - most_mean_diff_1/rowMeans(x_1a)
      df$fract_rmmd_d2gtd1[n] = sum(diff_rmmd > 0) / df$n_samples[n]
      df$mean_diff_rmmd_2m1[n] = mean(diff_rmmd)
      
      
      # standard normal random
      diff_nrand = rowMeans(matrix(rnorm(df$n_samples[n] * df$n_obs[n], mean = 0, sd = 1), 
                                   nrow = df$n_samples[n], ncol = df$n_obs[n]))
      df$fract_nrand_d2gtd1[n] = sum(diff_nrand > 0 ) / df$n_samples[n]
      df$mean_diff_nrand_2m1[n] = mean(diff_nrand)
      
      # Pearson Correlation
      # R Squared
      # Biserial Correlation https://rpubs.com/juanhklopper/biserial_correlation
    }
    # Save dataframed results to a file
    saveRDS(df, file = out_path)
  } else {
    # Restore the dataframed results from disk
    df <- readRDS(file = out_path)
  }
  return(df)
}



tidyEffectSizeResults <- function (df, gt_colname, var_suffix, long_format = TRUE,
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
      df_gt_sub[[matched_vars[n]]][ df[[gt_colname]]] <- 1-df[[matched_vars[n]]][ df[[gt_colname]]]
      df_gt_sub[[matched_vars[n]]][!df[[gt_colname]]] <-   df[[matched_vars[n]]][!df[[gt_colname]]]
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
      df_ref_sub[matched_vars[n]] <- df_gt_sub[matched_vars[n]]- df_gt_sub[ref_colname]
    }
  }
  
  # Flatten df_ref_sub into df_tidy. where each variable becomes a level of the factor 'Effect Size'
  if (long_format) {
    df_tidy <- df_ref_sub %>% gather(matched_vars,factor_key=TRUE)
    names(df_tidy)[names(df_tidy) == "matched_vars"] <- "name"
  } else df_tidy <- df_ref_sub
  
  return(df_tidy)
}

pretty_effect_size_levels<- function(df,base_names, pretty_names, var_suffix) {
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


# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomnly chosen
n_sims = 1e2
n_samples = 1e3
n_obs = 50
rand.seed = 0
fig_basename = "f_5"



# Contest 1) Quantify Error rate with each metric discerining experiment
# with lower mean difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = rep(1,n_sims), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_2a  = rep(1,n_sims), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_1b  = runif(n_sims, 1, 2), 
                                   sigmas_1b = rep(1,n_sims), 
                                   mus_2b  = runif(n_sims, 2, 3), 
                                   sigmas_2b = rep(1,n_sims)) 
# Quantify effect sizes in untidy matrix
df_mu <- quantifyEffectSizes(df = df_init,overwrite = TRUE, out_path = 
                               "temp/EffectSizeContest_mean_shift.rds")

# Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
df_mu_tidy <- tidyEffectSizeResults(df = df_mu, gt_colname = "is_sigmad_d2gtd1",
                                    var_suffix = var_suffix,long_format = TRUE,
                                    ref_colname = NULL)
df_mu_pretty <- pretty_effect_size_levels(df = df_mu_tidy, 
                                          base_names = paste(var_suffix,"_", effect_size_dict$base, "_", sep=""),
                                          pretty_names = effect_size_dict$label, 
                                          var_suffix = var_suffix)
# Basic violin plot
p <- ggplot(df_mu_pretty, aes(x=name,  y=value)) +
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  xlab("Effect Size Metric") +
  ylab(expression(Error~Rate~Exp.~Lower~mu[d])) +
  scale_x_discrete(labels = parse(text = levels(df_mu_pretty$name))) +
  theme_classic(base_size = 8) 
p
save_plot(paste("figure/", fig_basename, "1a effect size contest- mu.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.25, dpi = 600)


# Contest 2) Quantify Error rate with each metric discerining experiment
# with lower mean difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = runif(n_sims, -1, 1), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_2a  = runif(n_sims, -1, 1), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_1b  = runif(n_sims, 1, 2), 
                                   sigmas_1b = rep(1,n_sims), 
                                   mus_2b  = runif(n_sims, 1, 2), 
                                   sigmas_2b = rep(1,n_sims)) 
# Quantify effect sizes in untidy matrix
df_mu <- quantifyEffectSizes(df = df_init,overwrite = TRUE, out_path = 
                               "temp/EffectSizeContest_mean_shift.rds")

# Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
df_mu_tidy <- tidyEffectSizeResults(df = df_mu, gt_colname = "is_mud_d2gtd1",
                                    var_suffix = var_suffix,long_format = TRUE,
                                    ref_colname = NULL)#"fract_xdbar_d2gtd1")
df_mu_pretty <- pretty_effect_size_levels(df = df_mu_tidy, 
       base_names = paste(var_suffix,"_", effect_size_dict$base, "_", sep=""),
                                          pretty_names = effect_size_dict$label, 
                                          var_suffix = var_suffix)
# Basic violin plot
p <- ggplot(df_mu_pretty, aes(x=name,  y=value)) +
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  xlab("Effect Size Metric") +
  ylab(expression(Error~Rate~Exp.~Lower~mu[d])) +
  scale_x_discrete(labels = parse(text = levels(df_mu_pretty$name))) +
  theme_classic(base_size = 8) 
p
save_plot(paste("figure/", fig_basename, "1a effect size contest- mu.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.25, dpi = 600)




# Contest 2) which metric is better at discerining exp. with lower mean of
# relative difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_rmu_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                   mus_1a  = runif(n_sims, 10, 20), 
                                   sigmas_1a = rep(1,n_sims), 
                                   mus_2a  = runif(n_sims, 10, 20), 
                                   sigmas_2a = rep(1,n_sims),
                                   mus_1b  = runif(n_sims, 10, 20), 
                                   sigmas_1b = rep(1,n_sims), 
                                   mus_2b  = runif(n_sims, 10, 20), 
                                   sigmas_2b = rep(1,n_sims)) 
# Quantify effect sizes in untidy matrix
df_rmu <- quantifyEffectSizes(df = df_rmu_init,overwrite = TRUE, out_path = 
                               "temp/EffectSizeContest_mean_shift.rds")

# Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
df_rmu_tidy <- tidyEffectSizeResults(df = df_rmu, gt_colname = "is_rmud_d2gtd1",
                                    var_suffix = var_suffix,long_format = TRUE,
                                    ref_colname = NULL)#"fract_xdbar_d2gtd1")
df_rmu_pretty <- pretty_effect_size_levels(df = df_rmu_tidy, 
                                          base_names = paste(var_suffix,"_", effect_size_dict$base, "_", sep=""),
                                          pretty_names = effect_size_dict$label, 
                                          var_suffix = var_suffix)
# Basic violin plot
p <- ggplot(df_rmu_pretty, aes(x=name,  y=value)) +
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  xlab("Effect Size Metric") +
  ylab(expression(Error~Rate~Exp.~Lower~r*mu[d])) +
  scale_x_discrete(labels = parse(text = levels(df_mu_pretty$name))) +
  theme_classic(base_size = 8) 
p
save_plot(paste("figure/", fig_basename, "1a effect size contest- rmu.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.25, dpi = 600)




# Contest 3) which metric is better at discerining exp. with lower std of 
# difference in means
#
#------------------------------------------------------------------------------
var_suffix = "fract"
df_rmu_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                       mus_1a  = rep(0,n_sims), 
                                       sigmas_1a = rep(1,n_sims), 
                                       mus_2a  = rep(0,n_sims), 
                                       sigmas_2a = rep(1,n_sims),
                                       mus_1b  = rep(0,n_sims), 
                                       sigmas_1b = runif(n_sims, .1, 2),
                                       mus_2b  = rep(0,n_sims), 
                                       sigmas_2b = runif(n_sims, .1, 2))
# Quantify effect sizes in untidy matrix
df_rmu <- quantifyEffectSizes(df = df_rmu_init,overwrite = TRUE, out_path = 
                                "temp/EffectSizeContest_sigma.rds")

# Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
df_rmu_tidy <- tidyEffectSizeResults(df = df_rmu, gt_colname = "is_sigmad_d2gtd1",
                                     var_suffix = var_suffix,long_format = TRUE,
                                     ref_colname = NULL)#"fract_xdbar_d2gtd1")
df_rmu_pretty <- pretty_effect_size_levels(df = df_rmu_tidy, 
                                           base_names = paste(var_suffix,"_", effect_size_dict$base, "_", sep=""),
                                           pretty_names = effect_size_dict$label, 
                                           var_suffix = var_suffix)
# Basic violin plot
p <- ggplot(df_rmu_pretty, aes(x=name,  y=value)) +
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  xlab("Effect Size Metric") +
  ylab(expression(Error~Rate~Exp.~Lower~sigma[d])) +
  scale_x_discrete(labels = parse(text = levels(df_mu_pretty$name))) +
  theme_classic(base_size = 8) 
p
save_plot(paste("figure/", fig_basename, "1a effect size contest- sigma.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.25, dpi = 600)


# Contest 4) Quantify each metric's correlation with absolute relative effect size
#
#------------------------------------------------------------------------------




