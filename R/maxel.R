
# Load packages
#library(readr)
#library(dplyr)
library(ggplot2)
library(tibble)

nsamples = 15
ngroups = 2
set.seed(0)

# Sample two null distributions
null_dist = tibble(group = c(rep(0,nsamples), rep(1,nsamples)),
                   y = rnorm(ngroups*nsamples, 0, 1))




# Basic box plot
p <- ggplot(null_dist, aes(x=group,y=y)) + 
  geom_boxplot(aes(group=group,y=y)) + geom_jitter(width = .1)
# Fill in boxplot
p+scale_fill_brewer(palette="Dark2")
p




# 1B: distribution of effect between distributions

## 1A-C: One sample distribution case

## 1D-F: two sample distribution case

## 1G-I: 3 sample distirbution case











# MAXEL correlate with effect size

# MAXEL rewards experiments with lower noise

# MAXEL rewards experiment with lower bias

# MAXEL does not need a multiple test pvalue correction

