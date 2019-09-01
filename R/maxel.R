
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

equ_tbl = data.frame(maxel = c(expression(beta +frac(miles, gallon)),
                     expression(beta +frac(miles, gallon)),
                     expression(beta +frac(miles, gallon))), 
           rel_maxel = c(expression(beta +frac(miles, gallon)), 
                         expression(beta +frac(miles, gallon))),
           expression(beta +frac(miles, gallon)))


data.frame()
ggplot(as.data.frame(table(df)), aes(x=gender, y = Freq, fill=fraud)) + 
  geom_bar(stat="identity")



# Case 3: Unpaired data, two distributions
# 3.1 distribution A, B
x <- seq(-3,6,.01)
mu <- c(1,2)
data3_1 = rbind(
  tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=("2"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))

p3_1 <- ggplot(data3_1, aes(x=x, y=y, group=group, fill=group)) +
  geom_line(size=1,aes(color=group)) + 
  geom_ribbon(data3_1=subset(data3_1,x>-4 & x<7.5), 
              aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text", x=mu[1], y=0.4+.02, label= "mu[1]",
           parse=TRUE,color=color_pal[1],size=8) +
  annotate("text", x=mu[2], y=0.4+.02, label= "mu[2]",
           parse=TRUE,color=color_pal[2],size=8)

# 3.2: difference between distributions
x3_2 <- seq(-3,6,.01)
data3_2 = tibble(group=factor("1"), x=x3_2,
         y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
tt_3_2 = tidy(t.test(rnorm(1000, mean=mu[2]-mu[1], sd=1)))

p3_2 <-ggplot(data3_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data3_2$y,ymin=0, alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.02, label= "mu[d] == mu[a] - mu[b]",
           parse=TRUE,color=color_pal[3],size=8) +
  geom_vline(xintercept=tt_3_2$conf.low) +
  geom_vline(xintercept=tt_3_2$conf.high)



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

