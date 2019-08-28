
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



# Plot 2 distributions
p9 <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) +
  stat_function(fun = dt, args = list(df = 8))
p9


pdf_tbl = rbind(
  tibble(group="X1", x=seq(-5,5,.01) ,
         y=dnorm(seq(-5,5,.01), mean = 0, sd = 1, log = FALSE)),
  tibble(group="X2", x=seq(-5,5,.01) ,
         y=dnorm(seq(-5,5,.01), mean = 2, sd = 1, log = FALSE)))
p <- ggplot(pdf_tbl, aes(x=x,y=y,fill=group)) +
  geom_area(aes(fill=aes(group))) +
  geom_smooth(method="auto") 
p


p9 <- ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
  stat_function(fun = dnorm, args = list(0.2, 0.1),
                colour = color_pal[1]) +
  stat_function(fun = dnorm, args = list(0.7, 0.05),
                colour = color_pal[2]) +
  scale_x_continuous(name = "Probability",
                     breaks = seq(0, 1, 0.2),
                     limits=c(0, 1)) +
  scale_y_continuous(name = "Frequency") +
  ggtitle("Normal function curves of probabilities")
p9



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

