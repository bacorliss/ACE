

### Figure 2: MAXEL demonstrates a unique combination of properties compared to previous statistics of effect size.


### Measure of sample effect size
#  difference of sample means
#  cohen's d
#  g
#  Maxel

n_simulations=1e3;
n_trials = 1e6
n_observed = 10

### Maxel discerns experiments with lower variance
# Generate two matrices, X1 and X2 sampled from two standard normal distributions

# Generate 1,000 randon mean variances by coefficients sampled uniformly from -1:1

# Record difference in population variance

# Record whether sigma_1 > sigma_2 in binary vector

# Record whether MAXEL decreases in binoary vector

# Record whether cohen's d >1 

# Record whether cohen's g >1 

# Record whether sample variances decrease

# Plot accuracies calculated from each simulation 




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

