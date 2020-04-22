

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
library(gridExtra)
library(grid)
library(tidyr)
library(cowplot)
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
  function(m1,m2) sqrt( ( (dim(m1)[2] - 1) * rowSds(m1)^2   + 
                          (dim(m2)[2] - 1) * rowSds(m2)^2 ) /
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
  n1<-dim(m1)[2]
  n2<-dim(m2)[2] 
  n<-n1+n2 
  tstat<- sqrt(n1*n2/n) * ( rowMeans(m2) - rowMeans(m1)) /
          sqrt( ((n1 - 1) / (n - 2)* rowVars(m1) + (n2 - 1) / (n - 2)*rowVars(m2) )) 
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
effect_size_dict <- vector(mode="list", length=3)
names(effect_size_dict) <- c("prefix", "base", "suffix")
effect_size_dict[[1]] <- c("fract", "mean_diff")
effect_size_dict[[2]] <- c("means", "stds", "rel_means","rel_stds", "t_score",
                           "cohen_d", "hedge_g", "glass_delta", "most_mean_diff", 
                           "p_value")
effect_size_dict[[3]] <- c("2gt1","2m1")



generateExperiment_Data <- function(n_samples, n_obs, n_sims, rand.seed, 
                                    mus_1  = runif(n_sims, -0.5, 0.5), 
                                    sigmas_1 = runif(n_sims,  0.1, 1), 
                                    mus_2  = runif(n_sims, -0.5, 0.5), 
                                    sigmas_2 = runif(n_sims,  0.1, 1),
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

  # Output dataframe with altered
  set.seed(rand.seed +1)
  df = tibble(mu_1 = mus_1, sigma_1 = sigmas_1,
              mu_2 = mus_2, sigmas_2 = sigmas_2
  )
  # Is pop_mean2 > pop_mean1 (TRUE)
  df <- add_column(df, is_mu_d2gtd1 =  df$mu_2 > df$mu_1)
  #   pop_mean2 - pop_mean1
  df <- add_column(df, mean_mu_d2md1 = df$mu_2 - df$mu_1)
  # Is pop_std2 > pop_std1 (TRUE)
  df <- add_column(df, is_sigma_d2gtd1 =   df$sigma_2 > df$sigma_1)
  #   pop_std2 - pop_std1
  df <- add_column(df, mean_sigma_d2md1 =  df$sigma_2 - df$sigma_1)
  
  # Append columns for effect sizes, since multiple columns are used to analyze
  # each effect size, a dictionary of prefix,base, and suffix variable names 
  # are used.
  df[ paste(effect_size_dict$prefix[1], effect_size_dict$base, 
            effect_size_dict$suffix[1], sep="_") ] <- 0
  df[ paste(effect_size_dict$prefix[2], effect_size_dict$base, 
            effect_size_dict$suffix[2], sep="_") ] <- 0
  return(df)
}

quantifyEffectSizes <- function(df, x_a, x_b,out_path,overwrite=FALSE){
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
  #' 
  #' 
  
  # Only perform simulations if results not saved to disk
  if (!file.exists(out_path) | overwrite) {
    for (n in seq(1,dim(df)[1],1)) {
      # Transform simulated samples with normal parameteres (mean and std) for 
      # current round of simulation
      
      # Use Exp 1 and 2 coefficients to generate data from normalized base data
      x_1a = x_a
      x_1b = x_b * df$sigma_1[n] + df$mu_2[n]
      
      x_2a = x_a
      x_2b = x_b * df$sigma_2[n] + df$mu_1[n] 
      
      # Calculate difference of basic statistics
      dbar_1 = abs(rowMeans(x_1b) - rowMeans(x_1a))
      dbar_2 = abs(rowMeans(x_2b) - rowMeans(x_2a))
      
      sd_1 = rowSds(x_1b) - rowSds(x_1a)
      sd_2 = rowSds(x_2b) - rowSds(x_2b)
      
      # Basic Summary statistical comparisons
      # Means
      df$fract_dbar_2gtd1[n] = sum(dbar_2 > dbar_1)/n_samples
      df$mean_dbar_2md1[n]  = mean(dbar_2 - dbar_1)
      
      # Stds
      df$fract_sd_2gtd1[n] = sum(sd_2 > sd_1)/n_samples
      df$mean_diff_stds_2m1[n]  = mean(sd_2 - sd_1)
      
      # Rel Means: mean divided by control mean
      diff_rdbar = dbar_2/rowMeans(x_2a) - dbar_1/rowMeans(x_1a)
      df$fract_rel_means_2gt1[n] = sum( diff_rdbar > 0) / n_samples
      df$mean_diff_rel_means_2m1[n]  =  mean(diff_rdbar)
      
      
      # Rel STDs: sd divided by control mean
      diff_rsd = sd_2/rowMeans(x_2a) - sd_1/rowMeans(x_1a)
      df$fract_stds_2gt1[n] = sum( diff_rsd > 0) / n_samples
      df$mean_diff_stds_2m1[n]  = mean(diff_rsd)
      
      # Delta Family of effect size
      # Cohens D
      diff_cohen_d = abs(rowCohenD(x_2a, x_2b) - rowCohenD(x_1a, x_1b))
      df$fract_cohen_d_2gt1[n] = sum(diff_cohen_d >0) / n_samples
      df$mean_diff_cohen_d_2m1[n] =  mean(diff_cohen_d)
      # Glass delta
      diff_glass_delta = abs(rowGlassDelta(x_2a, x_2b) - rowGlassDelta(x_1a, x_1b))
      df$fract_glass_delta_2gt1[n] = sum(diff_glass_delta >0) / n_samples
      df$mean_diff_glass_delta_2m1[n] =  mean(diff_glass_delta)
      # Hedges G
      diff_hedge_g = abs(rowHedgeG(x_2a, x_2b) - rowHedgeG(x_1a, x_1b))
      df$fract_hedge_g_2gt1[n] = sum(diff_hedge_g >0) / n_samples
      df$mean_diff_hedge_g_2m1[n] =  mean(diff_hedge_g)
      
      # R Squared
      
      
      # Most Mean Diff
      most_mean_diff_1 = sapply(1:n_samples, function(i) mmd_normal(x_1a[i,], x_1b[i,]))
      most_mean_diff_2 = sapply(1:n_samples, function(i) mmd_normal(x_2a[i,], x_2b[i,]))
      diff_most_mean_diff = most_mean_diff_2 - most_mean_diff_1
      df$fract_most_mean_diff_2gt1[n] = sum(diff_most_mean_diff >0) / n_samples
      df$mean_diff_most_mean_diff_2m1[n] = mean(diff_most_mean_diff)
      
      # t score
      t_score_1 <- rowTScore(x_1a, x_1b)
      t_score_2 <- rowTScore(x_2a, x_2b)
      df$fract_t_score_2gt1[n] = sum( (t_score_2 - t_score_1) > 0) / n_samples
      df$mean_diff_most_mean_diff_2m1[n] = mean( (t_score_2 - t_score_1))
      
      # Pvalue
      diff_p_value <- pt(t_score_2, df = pmin(n_obs, n_obs) - 1) - 
        pt(t_score_1, df = pmin(n_obs, n_obs) - 1)
      df$fract_p_value_2gt1[n] = sum(diff_p_value > 0) / n_samples
      df$mean_diff_p_value_2m1 = mean(diff_p_value)
      
      # Pearson Correlation
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

tidyEffectSizeResults <- function (df, gt_vector, var_suffix) {
  # Get list of all variables that match the variable suffix string
  matched_vars <- grep(var_suffix,colnames(df), perl=TRUE, value=TRUE)
  
  # Initialize new df with only selected vars
  data = matrix(0, dim(df)[1], length(matched_vars))
  colnames(data) <- matched_vars
  df_select <- as_tibble(data)
 
  # Fill in resulting values depending on variable type
  # Frac: fraction of metrics of samples that satisfy some condition: to process,
  #        subtract these values from ground truth.
  # mean_diff: mean value of difference between groups:  to process,
  #        test for equaivalent logical values
  for (n in seq(1,length(matched_vars), by = 1)) {
    if (var_suffix == "fract") {
      # If var_suffix is "fract.*", for cases where GT vector is TRUE, compute (1-x)
      df_select[matched_vars[n]][gt_vector,] <- 1-df[matched_vars[n]][gt_vector,]
      df_select[matched_vars[n]][!gt_vector,] <- df[matched_vars[n]][!gt_vector,]
    } else {
      # if var_suffix is "mean_diff.*", subtract from GT
      df_select[matched_vars[n]] = df[matched_vars[n]] - gt_vector
    }
  }
  # Flatten df_select into df_tidy. where each variable becomes a level of the factor 'Effect Size'
  df_tidy <- df_select %>% gather(matched_vars)
}

# Simulation parameters
# 
#-------------------------------------------------------------------------------
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomnly chosen
n_sims = 1e2
n_samples = 1
n_obs = 50
rand.seed = 0

# Generate normalized measurement data
set.seed(rand.seed)
x_a = matrix(rnorm(n_samples*n_obs, mean = 0, sd = 1), n_samples, n_obs)
x_b = matrix(rnorm(n_samples*n_obs, mean = 0, sd = 1), n_samples, n_obs)


# Contest 1) Quantify Error rate with each metric discerining experiment
# with lower mean difference in means
#
#------------------------------------------------------------------------------
df_init <- generateExperiment_Data(n_samples, n_obs, n_sims, rand.seed, 
                                    mus_1  = rep(.1,n_sims), # runif(n_sims, -0.1, 0.1), 
                                    sigmas_1 = rep(1,n_sims), 
                                    mus_2  = rep(.1,n_sims), #runif(n_sims, -0.1, 0.1), 
                                    sigmas_2 = rep(1,n_sims)) 

# Quantify effect sizes in untidy matrix
df_mu_shifted <- quantifyEffectSizes(df_init, x_a, x_b,
                                     "temp/EffectSizeContest_mean_shift.rds", 
                                     overwrite=TRUE)
# Tidy matrix
df_mu_shifted_tidy <- tidyEffectSizeResults(df = df_mu_shifted, 
                                            gt_vector = df_mu_shifted$is_mu2d_gt_mu1d,
                                            var_suffix = "fract")
# Basic violin plot
p <- ggplot(df_mu_shifted_tidy, aes(x=matched_vars,  y=value)) +
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9)) +
  theme_classic(base_size = 8) + theme(legend.position="none",
                                       axis.title.x = element_blank()) +
  xlab("Metric") +
  ylab("Error Rate Exp with Lower mu")
p

save_plot(paste("figure/", fig_basename, "f mmd transition values.tiff",
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)





# Contest 2) which metric is better at discerining exp. with lower std of 
# difference in means
#
#------------------------------------------------------------------------------





# Contest 3) which metric is correlates best with absolute effect size
#
#------------------------------------------------------------------------------



# Contest 4) Quantify each metric's correlation with absolute relative effect size
#
#------------------------------------------------------------------------------




