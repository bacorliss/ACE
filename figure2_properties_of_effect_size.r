

### Figure 2: MAXEL demonstrates a unique combination of properties compared to previous statistics of effect size.



library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)

library(lsr)

# User defined libraries
source("R/maxel.R")

# USer defined functions
# Calculate variance by row like rowMeans or rowSums
# sum(x-x_bar).^2/(n-1)
rowVars <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)}
cohensD <- function(x,n,...) {
  (rowMeans(x, (n + 1):(2 * n)) - rowMeans(x, 1 : n))
}


### Measure of sample effect size
#  Unstandardized measures of effect size
#  difference of sample means
# Standardized measures of effect size
#  cohen's d
#  g
#  Maxel

n_sims=1e3;
n_samples = 1e4
n_obs = 10

# A simulation is a series of trials, each with 2 groups of 10 observations

### Maxel discerns experiments with lower variance
# Generate two matrices, X1 and X2 sampled from two standard normal distributions
set.seed(0)
exp_1 = matrix(rnorm(2*n_samples*n_obs, mean=1,sd=1),n_samples,2*n_obs)
exp_2 = matrix(rnorm(2*n_samples*n_obs, mean=1,sd=1),n_samples,2*n_obs)
sim_means_coeff = runif(n_sims,-0.5,0.5)

#sweep(b, 2, a, "+"

# Output dataframe with altered
df = tibble(pop_mean1 = rep(1,n_sims), pop_mean2 = 1 + sim_means_coeff,
                    pop_var1=rep(1,n_sims),pop_var2=rep(1,n_sims));
# Groundtruth: is sigma_2 > sigma_1 (= 1)
df <- add_column(df, is2gt1_pop_mean= df$pop_mean1 > df$pop_mean2)
# Groundtruth: difference of sigma2 - sigma_1 (= 1)
df <- add_column(df, diff2m1_pop_mean= df$pop_mean2 - df$pop_mean2)
# For each effect size, need 
#     frac_2gt1_[es]
#     mean_2m1_[es]

# Get list of effect size types

# Append effect size suffixes to all effect sizes

# Add all fields to refults df initialized to 0
df[c("mean_2m1_sample_mean", "frac_2gt1_sample_mean",
      "mean_2m1_sample_var", "frac_2gt1_sample_var")] <- 0


# Loop through simulations
# for (n in seq(1,n_sims,1)) {}

n=1
# Multiplly half of orig_x1_2 by pop_var2
x1_2 = sweep(orig_x1_2, MARGIN=2, c(rep(0, n_obs), rep(df$pop_var2[n], n_obs)), `*`)

# Calculate mean difference of sample means
sample_means = cbind(rowMeans(x1_2[ ,1:n_obs]), rowMeans(x1_2[ ,(n_obs+1) : (2*n_obs)]))
# Calculate mean value and fraction of group 2 > group 1
df$mean_2m1_sample_mean[n] = mean(sample_means[,2]-sample_means[,1])
df$frac_2gt1_sample_mean[n] = sum(sample_means[,2] > sample_means[,1])/n_samples


# Calculate 


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

