

### Figure 2: MHD outperforms previous statistics of effect size.



library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)
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
# Effect Size Statistics
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

# Simulation parameters
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomnly chosen
n_sims = 1e2
n_samples = 1e4
n_obs = 50
rand.seed = 0

### Maxel discerns experiments with lower variance
# Generate two matrices, X1 and X2 sampled from two standard NORMalized distributions
#  X1 is designated as the control sample that is used to normalize values for 
#  relative statistical metrics


# Code structure: 
# 1. two groups of untransformed sample sets are generated (x_a,x_b)
# 2. Transform coefficients for experiment 1 and 2 are generated
# 3. For each simulation, coefficients are applied
#     1. x_a [exp_1_coeffs]-> x1a
#     2. x_b [exp_1_coeffs]-> x1b
#     1. x_a [exp_2_coeffs]-> x2a
#     2. x_b [exp_2_coeffs]-> x2b
#     
#                     x_a                  x_b
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

set.seed(rand.seed)
x_a = matrix(rnorm(n_samples*n_obs, mean = 0, sd = 1), n_samples, n_obs)
x_b = matrix(rnorm(n_samples*n_obs, mean = 0, sd = 1), n_samples, n_obs)


# Output dataframe with altered
set.seed(rand.seed +1)
df = tibble(exp1_pop_mean = runif(n_sims, -0.5, 0.5), 
            exp1_pop_std = runif(n_sims, 0.1, 1),
            exp2_pop_mean = runif(n_sims, -0.5, 0.5), 
            exp2_pop_std = runif(n_sims, 0.1, 1))
# Is pop_mean2 > pop_mean1 (TRUE)
df <- add_column(df, pop_mean_is_2gt1 = df$exp2_pop_mean > df$exp1_pop_mean)
#   pop_mean2 - pop_mean1
df <- add_column(df, pop_mean_diff_2m1 = df$exp2_pop_mean - df$exp1_pop_mean)
# Is pop_std2 > pop_std1 (TRUE)
df <- add_column(df, pop_std_is_2gt1 =   df$exp2_pop_std > df$exp1_pop_std)
#   pop_std2 - pop_std1
df <- add_column(df, pop_std_diff_2m1 = df$exp2_pop_std - df$exp1_pop_std)

# Generate effect size columns in datatable
#   For each effect size metric, need 
#     frac_2gt1_[metric]
#     mean_2m1_[metric]
effect_size_prefix = c("fract", "mean_diff")
effect_sizes = c("means", "stds", "rel_means","rel_stds", "t_score",
                    "cohen_d", "hedge_g", "glass_delta", "most_mean_diff", "p_value")
# Append effect size suffixes to all effect sizes
df[ paste(effect_size_prefix[1], effect_sizes, "2gt1", sep="_") ] <- 0
df[ paste(effect_size_prefix[2], effect_sizes, "2m1", sep="_") ] <- 0


# Loop through simulations, using different set of parameters
for (n in seq(1,n_sims,1)) {

  # Transform simulated samples with normal parameteres (mean and std) for 
  # current round of simulation
  
  # Use Exp 1 and 2 coefficients to generate data
  x_1a = x_a * df$exp1_pop_std[n] + df$exp1_pop_mean[n]
  x_1b = x_b * df$exp1_pop_std[n] + df$exp1_pop_mean[n]
  x_2a = x_a * df$exp2_pop_std[n] + df$exp2_pop_mean[n] 
  x_2b = x_b * df$exp2_pop_std[n] + df$exp2_pop_mean[n] 
  
  # Calculate difference of basic statistics
  exp1_diff_means = rowMeans(x_1b) - rowMeans(x_1a)
  exp2_diff_means = rowMeans(x_2b) - rowMeans(x_2a)
  exp1_diff_sds = rowSds(x_1b) -rowSds(x_1a)
  exp2_diff_sds = rowSds(x_2b) -rowSds(x_2b)

  # Basic Summary statistical comparisons
  # Means
  df$fract_means_2gt1[n] = sum(exp2_diff_means > exp1_diff_means)/n_samples
  df$mean_diff_means_2m1[n]  = mean(exp2_diff_means - exp1_diff_means)
  # Stds
  df$fract_stds_2gt1[n] = sum(exp2_diff_sds > exp1_diff_sds)/n_samples
  df$mean_diff_stds_2m1[n]  = mean(exp2_diff_sds - exp1_diff_sds)
  # Rel Means: mean divided by control mean
  diff_rel_mean_diff = exp2_diff_means/rowMeans(x_2a) - exp1_diff_means/rowMeans(x_1a)
  df$fract_rel_means_2gt1[n] = sum( diff_rel_mean_diff > 0) / n_samples
  df$mean_diff_rel_means_2m1[n]  =  mean(diff_rel_mean_diff)
  # Rel STDs: sd divided by control mean
  diff_rel_sds_diff = exp2_diff_sds/rowMeans(x_2a) - exp1_diff_sds/rowMeans(x_1a)
  df$fract_stds_2gt1[n] = sum( diff_rel_sds_diff > 0) / n_samples
  df$mean_diff_stds_2m1[n]  = mean(diff_rel_sds_diff)

  # Delta Family of effect size
  # Cohens D
  diff_cohen_d = rowCohenD(x_2a, x_2b) - rowCohenD(x_1a, x_1b)
  df$fract_cohen_d_2gt1[n] = sum(diff_cohen_d >0) / n_samples
  df$mean_diff_cohen_d_2m1[n] =  mean(diff_cohen_d)
  # Glass delta
  diff_glass_delta = rowGlassDelta(x_2a, x_2b) - rowGlassDelta(x_1a, x_1b)
  df$fract_glass_delta_2gt1[n] = sum(diff_glass_delta >0) / n_samples
  df$mean_diff_glass_delta_2m1[n] =  mean(diff_glass_delta)
  # Hedges G
  diff_hedge_g = rowHedgeG(x_2a, x_2b) - rowHedgeG(x_1a, x_1b)
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
### Maxel discerns results with lower bias




### Maxel reflects the true effect size between population means

# Perform 100,000 t-test for 1,000 different effect sizes

# Generate two matrices, X1 and X2 sampled from two standard normal distributions
## rows are observations, and columns are sampled trials
# Gernate 1,000 random mean offsets from a uniform distribution -2:2


# Calculate population effect size between populations, u1 and u2

# For each trial, calculate various measure of sample effect size


# Perform regression/ (multiple regression?) to assess predictive power for each effect size


# Visualize accuracy/performance for each model
# Show that MAXEL performs nearly as well as difference of sample means

















# plot(1,1, main=expression('title'[2]))

#=max(abs(ub-xbar, xbar-lb))

## 1B: distribution of effect between distributions

## 1A-C: One sample distribution case

## 1D-F: two sample distribution case

## 1G-I: 3 sample distirbution case


# Sample two null distributions
null_dist = tibble(group = factor(c(rep("x1",nsamples), rep("x2",nsamples))),
                   y = rnorm(ngroups*nsamples, 0, 1)+1)

# Basic box plot
p <- ggplot(null_dist, aes(x=group,y=y,fill=group)) + 
  geom_boxplot(aes(group=group,y=y),notch=FALSE) + 
  stat_boxplot(geom = 'errorbar')+
  geom_jitter(width = .1) +
  scale_fill_manual(values = color_pal[1:2]) + 
  scale_x_discrete(labels = c(expression(X[1]),expression(X[2]))) +
  theme(legend.position = "none")
p
# calculate ttest result
t_test_result <- t.test(y~group, null_dist,
                        alternative = "two.sided",
                        mu = 0,
                        paired = FALSE,
                        var.equal = FALSE,
                        conf.level = 0.95) %>% tidy()








# MAXEL correlate with effect size

# MAXEL rewards experiments with lower noise

# MAXEL rewards experiment with lower bias

# MAXEL does not need a multiple test pvalue correction

# MAXEL's value can be used for interpreting equivalence

