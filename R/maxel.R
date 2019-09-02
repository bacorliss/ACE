
# Load packages
#library(readr)
#library(dplyr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)

# choose colors for plotting
color_pal = brewer.pal(3, "Set1")

# Initial parameters
nsamples = 15
ngroups = 2
set.seed(0)
dist_text_size=3;


### Case 1: one-sample
# 1.1 distribution A, B
x <- seq(-3,6,.01)
mu <- c(1,2)
sd <- c(1,1)
data1_1 = rbind(
  tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=("2"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))

p1_1 <- ggplot(data1_1, aes(x=x, y=y, group=group, fill=group)) +
  geom_line(size=1,aes(color=group)) + 
  geom_ribbon(data=subset(data1_1,x>-4 & x<7.5), 
              aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text", x=mu[1], y=0.4+.02, label= "mu[1]",
           parse=TRUE,color=color_pal[1],size=dist_text_size) +
  annotate("text", x=mu[2], y=0.4+.02, label= "mu[2]",
           parse=TRUE,color=color_pal[2],size=dist_text_size)

### 3.2: difference with one-sample
x3_2 <- seq(-3,6,.01)
data1_2 = tibble(group=factor("1"), x=x3_2,
                 y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
p1_2 <-ggplot(data, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data1_2$y,ymin=0, alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.03, label= 
             "d == N(mu[A]-mu[B],sigma^2[A] + sigma^2[B])",
           parse=TRUE,color=color_pal[3],size=dist_text_size) +
  geom_vline(xintercept=mu[2]-mu[1]-2*sd[1]) +
  geom_vline(xintercept=mu[2]-mu[1]+2*sd[2])



### Case 2: two sample paired
# 2.1 distribution A, B
x <- seq(-3,6,.01)
mu <- c(1,2)
sd <- c(1,1)
data2_1 = rbind(
  tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=("2"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))

p2_1 <- ggplot(data2_1, aes(x=x, y=y, group=group, fill=group)) +
  geom_line(size=1,aes(color=group)) + 
  geom_ribbon(data=subset(data2_1,x>-4 & x<7.5), 
              aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text", x=mu[1], y=0.4+.02, label= "mu[1]",
           parse=TRUE,color=color_pal[1],size=dist_text_size) +
  annotate("text", x=mu[2], y=0.4+.02, label= "mu[2]",
           parse=TRUE,color=color_pal[2],size=dist_text_size)

### 3.2: difference two sample paired
x3_2 <- seq(-3,6,.01)
data2_2 = tibble(group=factor("1"), x=x3_2,
                 y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
p2_2 <-ggplot(data2_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data2_2$y,ymin=0, alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.03, label= 
             "d == N(mu[A]-mu[B],sigma^2[A] + sigma^2[B])",
           parse=TRUE,color=color_pal[3],size=dist_text_size) +
  geom_vline(xintercept=mu[2]-mu[1]-2*sd[1]) +
  geom_vline(xintercept=mu[2]-mu[1]+2*sd[2])


### Case 3: Unpaired two-sample
# 3.1 distribution A, B
x <- seq(-3,6,.01)
mu <- c(1,2)
sd <- c(1,1)
data3_1 = rbind(
  tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=("2"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))

p3_1 <- ggplot(data3_1, aes(x=x, y=y, group=group, fill=group)) +
  geom_line(size=1,aes(color=group)) + 
  geom_ribbon(data=subset(data3_1,x>-4 & x<7.5), 
              aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text", x=mu[1], y=0.4+.02, label= "mu[1]",
           parse=TRUE,color=color_pal[1],size=dist_text_size) +
  annotate("text", x=mu[2], y=0.4+.02, label= "mu[2]",
           parse=TRUE,color=color_pal[2],size=dist_text_size)

### 3.2: difference of Unpaired two-sample
x3_2 <- seq(-3,6,.01)
data3_2 = tibble(group=factor("1"), x=x3_2,
         y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
p3_2 <-ggplot(data3_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data3_2$y,ymin=0, alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.03, label= 
             "d == N(mu[A]-mu[B],sigma^2[A] + sigma^2[B])",
           parse=TRUE,color=color_pal[3],size=dist_text_size) +
  geom_vline(xintercept=mu[2]-mu[1]-2*sd[1]) +
  geom_vline(xintercept=mu[2]-mu[1]+2*sd[2])


### Figure creation
## Assign objects to GROB table
# Column names
gs[[1]] <- textGrob("Distribution", x = unit(0.45, "npc"), 
                    y = unit(0.5, "npc"), just = "centre",
                    gp=gpar(fontface="bold"))
gs[[2]] <- textGrob("Effect", x = unit(0.5, "npc"), 
                    y = unit(0.5, "npc"), just = "centre",
                    gp=gpar(fontface="bold"))
gs[[3]] <- textGrob("MAXEL", x = unit(0.5, "npc"), 
                    y = unit(0.5, "npc"), just = "centre",
                    gp=gpar(fontface="bold"))
gs[[4]] <- textGrob("Rel. MAXEL", x = unit(0.5, "npc"), 
                    y = unit(0.5, "npc"), just = "centre",
                    gp=gpar(fontface="bold"))

# One sample test
gs[[5]] <- p1_1
gs[[6]] <- p1_2
gs[[7]] <- textGrob("d == bar(X[A]) + bar(X[B])", x = unit(0.5, "npc"), 
                    y = unit(0.5, "npc"), just = "centre")
gs[[8]] <- textGrob("d == bar(X[A]) + bar(X[B])", x = unit(0.5, "npc"), 
                       y = unit(0.5, "npc"), just = "centre")
# Two sample test, paired
gs[[9]] <- p2_1
gs[[10]] <- p2_2
gs[[11]] <- textGrob("d == bar(X[A]) + bar(X[B])", x = unit(0.5, "npc"), 
                        y = unit(0.5, "npc"), just = "centre")
gs[[12]] <- textGrob("d == bar(X[A]) + bar(X[B])", x = unit(0.5, "npc"), 
                     y = unit(0.5, "npc"), just = "centre")
# Two sample test, unpaired
gs[[13]] <- p3_1
gs[[14]] <- p3_2
gs[[15]] <- textGrob("d == bar(X[A]) + bar(X[B])", x = unit(0.5, "npc"), 
                     y = unit(0.5, "npc"), just = "centre")
gs[[16]] <- textGrob("d == bar(X[A]) + bar(X[B])", x = unit(0.5, "npc"), 
                     y = unit(0.5, "npc"), just = "centre")
# Arrange grob obejcts into grid
gt <- arrangeGrob(grobs = gs, layout_matrix = rbind(c(1,2,3,4),
                                               c(5,6,7,8),
                                               c(9,10,11,12),
                                               c(13,14,15,16)))

gt <- gtable::gtable_add_grob(gt, 
                             grobs = rectGrob(gp=gpar(fill=NA, 
                                                      lwd=2)), 
                             t = 1, b = 1, l = 1, r = ncol(gt))

# Change relative height of gtable rows
gt$heights <- unit(c(.05, .3, .3, .3), "npc")
# Render new figure
grid.newpage()
grid.draw(gt)








# Table of MAXEL equations
equ_tbl = data.frame(maxel = c("d == bar(X[A]) + bar(X[B])",
                               "d == bar(X[A]) + bar(X[B])",
                               "d == bar(X[A]) + bar(X[B])"), 
                     rel_maxel = c("d==bar(X[A])+bar(X[B])", 
                                   "d == bar(X[A]) + bar(X[B])",
                                   "d == bar(X[A]) + bar(X[B])"))
colnames(equ_tbl) <- c("MAXEL", "Rel. MAXEL")
tt <- ttheme_default(core=list(fg_params = list(parse=TRUE)))
g <- tableGrob(equ_tbl, theme=tt)
g$heights <- unit(rep(1/nrow(g), nrow(g)), "npc")
grid.draw(g)
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

