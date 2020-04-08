

### Figure 2: MHD outperforms previous statistics of effect size.



library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)
library(effsize)
library(lsr)

# User defined libraries
source("R/maxel.R")

## User defined functions
# Calculate variance by row like rowMeans or rowSums
rowVars <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)}
rowSds  <- function(x, ...) sqrt(rowVars(x))
# Calculate row-wise Cohens D
sdPooled <- function(x1,x2) sqrt( ( (length(x1)-1) * sd(x1)^2 + 
                                    (length(x2)-1) * sd(x2)^2 ) /
                                    (length(x1) + length(x2) - 2))
# Effect Size Statistics
#   d = (M2 - M1)/s_pool      s_pool = sqrt( (s1 + s2)/n)
cohensD <- function (x1,x2) cohen.d(x1,x2,pooled=TRUE,paired=FALSE,
                                    na.rm=FALSE, mu=0, hedges.correction=FALSE,
                                    conf.level=0.95,noncentral=FALSE, ...)
#   d = (M2 - M1)/s_pool      s_pool = sqrt( (s1 + s2)/n)
hedgesG <- function (x1,x2) cohen.d(x1,x2,pooled=TRUE,paired=FALSE,
                                    na.rm=FALSE, mu=0, hedges.correction=TRUE,
                                    conf.level=0.95,noncentral=FALSE, ...)
#   d = (M2 - M1)/s_1
glassDelta <- function (x1,x2) (mean(x2)- mean(x1))/sd(x1) 

# 
# rowCohensD <- function(x,n,...) {
#   (rowMeans(x, (n + 1):(2 * n)) - rowMeans(x, 1 : n))
# }


# Metrics of sample effect size
##############################################
#  Unstandardized measures of effect size
#    difference of sample means
#    relative difference of means
#    difference of sample variances
#
#  Standardized measures of effect size
#    cohen's d
#    hedges g
#    glasses delta
#    coefficient of variation
#    MHD

# Simulation parameters
# A simulation is a set of samples with a fixed set of parameters
# Parameters are randomnly chosen
n_sims = 1e3;
n_samples = 1e4
n_obs = 10


### Maxel discerns experiments with lower variance
# Generate two matrices, X1 and X2 sampled from two standard NORMalized distributions
#  X1 is designated as the control sample that is used to normalize values for 
#  relative statistical metrics
set.seed(0)
norm_x1 = matrix(rnorm(n_samples*n_obs, mean = 0, sd = 1), n_samples, n_obs)
norm_x2 = matrix(rnorm(n_samples*n_obs, mean = 0, sd = 1), n_samples, n_obs)


# Output dataframe with altered
df = tibble(pop_mean1 = runif(n_sims, -0.5, 0.5), 
            pop_mean2 = runif(n_sims, -0.5, 0.5),
            pop_std1 = runif(n_sims, 0.01, 1), 
            pop_std2 = runif(n_sims, 0.01, 1));
# Is pop_mean2 > pop_mean1 (TRUE)
df <- add_column(df, pop_mean_is_2gt1 = df$pop_mean2 > df$pop_mean1)
#   pop_mean2 - pop_mean1
df <- add_column(df, pop_mean_diff_2m1 = df$pop_mean2 - df$pop_mean1)
# Is pop_std2 > pop_std1 (TRUE)
df <- add_column(df, pop_std_is_2gt1 =   df$pop_std2 > df$pop_std1)
#   pop_std2 - pop_std1
df <- add_column(df, pop_std_diff_2m1 = df$pop_std2 - df$pop_std1)


# Generate effect size columns in datatable
#   For each effect size metric, need 
#     frac_2gt1_[metric]
#     mean_2m1_[metric]
effect_size_prefix = c("fract_2g1", "diff_2m1")
effect_sizes = c("means", "stds", "rel_means","rel_stds",
                    "cohen_d", "hedge_g", "glass_delta", "most_mean_diff", "p_value", "coeff_var")
# Append effect size suffixes to all effect sizes
df[paste(effect_size_prefix[1], effect_sizes, sep="_")] <- 0
df[paste(effect_size_prefix[2], effect_sizes, sep="_")] <- 0


# Loop through simulations, using different set of parameters
for (n in seq(1,n_sims,1)) {

  # Transform simulated samples with normal parameteres (mean and std) for 
  # current round of simulation
  
  # Multiplly half of orig_x1_2 by pop_var2
  x1 = norm_x1 * df$pop_std1[n] + df$pop_mean1[n]
  x2 = norm_x2 * df$pop_std2[n] + df$pop_mean2[n] 
  means_x1 = rowMeans(x1)
  means_x2 = rowMeans(x2)
  sds_x1 = rowSds(x1)
  sds_x2 = rowSds(x2)
  
  # Calculate mean difference of sample means
  sample_means = cbind(rowMeans(x1_2[ ,1:n_obs]), rowMeans(x1_2[ ,(n_obs+1) : (2*n_obs)]))
  # Basic Summary statistical comparisons
  # Means
  df$frac_means_2gt1[n] = sum(means_x2 > means_x1)/n_samples
  df$diff_means_2m1[n]  = mean(means_x2 - means_x1)
  # Stds
  df$frac_stds_2gt1[n] = sds_x2 > sds_x1
  df$diff_stds_2m1[n]  = sds_x2 - sds_x1
  # Rel Means
  df$frac1_rel_means_2gt1[n] = sum( (means_x2 - means_x1)/means_x1 > 0)/n_samples
  df$diff_rel_means_2m1[n]  = mean((means_x2 - means_x1)/means_x1)
  # Rel STDs
  df$frac_stds_2gt1[n] = sds_x2 > sds_x1
  df$diff_stds_2m1[n]  = sds_x2 - sds_x1
  # Cohens D
  df$frac_cohen_d_2gt1[n] = cohensD(x1, x2) > 0
  df$cohen_d_2m1[n] =  cohensD(x1, x2)
  # Glass delta
  df$frac_glass_delta_2gt1[n] = cohensD(x1, x2) > 0
  df$glass_delta_2m1_[n] =  cohensD(x1, x2)
  # Hedges G
  df$frac_hedges_g_2gt1[n] = cohensD(x1, x2) > 0
  df$hedges_g_2m1 =  cohensD(x1, x2)
  # Most Mean Diff
  df$frac_most_mean_diff_2gt1[n] = cohensD(x1, x2) > 0
  df$most_mean_diff_2m1[n] = 
    # Pvalue
  df$frac_p_value_2gt1[n] = t.test(x1, x2) > 0.05
  df$p_value_2m1 = hedgesG(x1, x2)
  # Coefficient of variation
  df$frac_coeff_var_2gt1[n] = sds_x2/mean_x2 > sds_x1/mean_x1
  df$coeff_var_2m1[n] = sds_x2/mean_x2 - sds_x1/mean_x1
  
}
### Maxel discerns results with lower bias


df$mean_maxel = mean(apply(x1_2, 2, function(x) 
  maxel.normal_unpairedt(x=x1_2[1:n_obs],y=x1_2[(n_obs+1):(2*n_obs)]) ))

# Compute effect sizes on a row by row basis
temp_df = compute_effect_size(x = x1_2[, 1:n_obs],y = x1_2[(n_obs+1):(2*n_obs)]) 
# Extract measures of effect size
df$mean_diff_of_means[n] = mean(temp_df$difference_of_means)
df$frac_gt0_diff_of_means[n] = sum(temp_df$difference_of_means > 0)/ n_samples
df$mean_diff_of_variances[n] = mean(temp_df$difference_of_variances)
df$frac_gt0_diff_of_variances[n] = sum(temp_df$difference_of_variances > 0)/ n_samples
df$mean_cohens_d[n] = mean(temp_df$cohens_d)
df$frac_gt0_cohens_d[n] = sum(temp_df$cohens_d > 0)/ n_samples
df$mean_glass_delta[n] = mean(temp_df$glass_delta) 
df$frac_gt0_glass_delta[n] = sum(temp_df$glass_delta > 0)/ n_samples
df$mean_hedges_g[n] = mean(temp_df$hedges_g)
df$frac_gt0_hedges_g[n] = sum(temp_df$hedges_g > 0)/ n_samples


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

