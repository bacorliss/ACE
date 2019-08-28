
# Load packages
#library(readr)
#library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)

color_pal = brewer.pal(3, "Set1")

nsamples = 15
ngroups = 2
set.seed(0)



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

# plot(1,1, main=expression('title'[2]))

#=max(abs(ub-xbar, xbar-lb))

## 1B: distribution of effect between distributions

## 1A-C: One sample distribution case

## 1D-F: two sample distribution case

## 1G-I: 3 sample distirbution case










# MAXEL correlate with effect size

# MAXEL rewards experiments with lower noise

# MAXEL rewards experiment with lower bias

# MAXEL does not need a multiple test pvalue correction

