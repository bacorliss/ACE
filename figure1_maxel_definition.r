
# This script produces a figure defining MAXEL for one-sample and two-sample two-tailed datasets
# Bruce Corliss
# 9/1/2019

# Load packages
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)

library(rlang)

setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

# choose colors for plotting
color_pal = brewer.pal(3, "Set1")

# Initial parameters
nsamples = 15
ngroups = 2
set.seed(0)
dist_text_size=5;
title_font_size=7;
equ_font_size=8;

### Case 1: one-sample
# 1.1 distribution A, B
x <- seq(-3,6,.01)
mu <- 2
sd <- 1
data1_1 = tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu, sd = sd, log = FALSE))$y)
C0 <- 1;
p1_1 <- ggplot(data1_1, aes(x = x, y = y, group = group, fill = color_pal[2])) +
  geom_line(size = 1, color = color_pal[2]) + 
  geom_ribbon(data = data1_1, aes(x = x,ymax = y),ymin = 0, alpha=0.3, color = color_pal[2], fill = color_pal[2]) +
  annotate("text", x = 5, y = 0.4 + 0.02, label = "A",
           parse=TRUE,color=color_pal[2],size = dist_text_size) + 
  geom_vline(xintercept=C0, color = color_pal[1]) + 
  annotate("text", x = 5, y = 0.4 + 0.02, label = "C[0]",
           parse=TRUE,color=color_pal[1],size = dist_text_size) +
  theme(plot.title = element_text(hjust = 1,size = title_font_size),legend.position = "none") 


### 1.2: difference with one-sample
data2_1 <- data1_1
data2_1$x <- data1_1$x-C0
p1_2 <- ggplot(data2_1, aes(x = x, y = y, group = group, fill = color_pal[3])) +
  geom_line(size = 1, color = color_pal[3]) + 
  geom_ribbon(data = data2_1, aes(x = x,ymax = y),ymin = 0, alpha=0.3, 
              color = color_pal[2], fill = color_pal[3]) +
  annotate("text", x = 4, y = 0.4 + 0.02, label = "d",
           parse=TRUE,color=color_pal[3],size = dist_text_size) + 
#  geom_vline(xintercept=C0, color = color_pal[1]) + 
#  annotate("text", x = C0 - .3, y = 0.4 + 0.02, label = "C[0]",
#           parse=TRUE,color=color_pal[1],size = dist_text_size) +
  theme(plot.title = element_text(hjust = 1, size = title_font_size)) +
  geom_vline(xintercept=-2*sd) +
  geom_vline(xintercept=+2*sd)


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
  geom_ribbon(data=data2_1, 
              aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text", x=mu[1], y=0.4+.02, label= "mu[1]",
           parse=TRUE,color=color_pal[1],size=dist_text_size) +
  annotate("text", x=mu[2], y=0.4+.02, label= "mu[2]",
           parse=TRUE,color=color_pal[2],size=dist_text_size) +
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none") 


### 2.2: difference two sample paired
x3_2 <- seq(-3,6,.01)
data2_2 = tibble(group=factor("1"), x=x3_2,
                 y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
p2_2 <-ggplot(data2_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data2_2$y,aes(ymin=0), alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.03, label = "d",
           parse=TRUE,color=color_pal[3],size=dist_text_size) +
  geom_vline(xintercept=mu[2]-mu[1]-2*sd[1]) +
  geom_vline(xintercept=mu[2]-mu[1]+2*sd[2]) +
  theme(plot.title = element_text(hjust = 1, size = title_font_size)) 

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
  geom_ribbon(data=data3_1, 
              aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text", x=mu[1], y=0.4+.02, label= "mu[1]",
           parse=TRUE,color=color_pal[1],size=dist_text_size) +
  annotate("text", x=mu[2], y=0.4+.02, label= "mu[2]",
           parse=TRUE,color=color_pal[2],size=dist_text_size) + 
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none") 


### 3.2: difference of Unpaired two-sample
x3_2 <- seq(-3,6,.01)
data3_2 = tibble(group=factor("1"), x=x3_2,
                 y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
p3_2 <-ggplot(data3_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data3_2$y,aes(ymin=0), alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.03, label= 
             "d",parse=TRUE,color=color_pal[3],size=dist_text_size) +
  geom_vline(xintercept=mu[2]-mu[1]-2*sd[1]) +
  geom_vline(xintercept=mu[2]-mu[1]+2*sd[2]) +
  theme(plot.title = element_text(hjust = 1, size = 10)) 


### Figure creation
## Assign objects to GROB table
# Column names
gs <- lapply(1:16, function(ii) grobTree(rectGrob(gp = gpar(fill = ii, alpha = 0.5)), textGrob(ii)))
gs[[2]] <- textGrob("Initial Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[3]] <- textGrob("Effect Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[4]] <- textGrob("Equations", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))


# Row Names for figure
gs[[1]] <- textGrob("", rot = 90)
gs[[5]] <- textGrob("1-Sample", rot = 90, just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[9]] <- textGrob("2-Sample, Paired", rot = 90, just = "centre", gp=gpar(fontface = "bold", fontsize = 10))
gs[[13]] <- textGrob("2-Sample, Unpaired", rot = 90, just = "centre", gp = gpar(fontface = "bold", fontsize = 10))

### Equations for figure %>% 


df_1s <- data.frame(equations = c(
  "A==N(mu[A]==2,sigma[A]==1)",
  "C[0]==1",
  "italic(d) == N(mu[A]-C[0],{sigma^2}[A])",
  "MXL==max*({}~abs({}~italic(bar(x)[1])-beta%+-%1.96~s[1]~{})~{})"
  ))

df_2s_p <- data.frame(equations = c(
  "A==N(mu[A]==1,sigma[A]==1)",
  "B==N(mu[B]==2,sigma[B]==1)",
  "d[i]==x[i,a]-x[i,b]",
  "MXL==max*({}~abs({}~italic(bar(d))%+-%1.96~s[diff]~{})~{})",
  "Rel~MXL==frac(max*({}~abs({}~italic(bar(d))%+-%1.96~s[diff]~{})~{}),min*({}~abs({}~italic(bar(d))%+-%italic(t)[list(0.975,n-1)]~frac(s[diff],sqrt(n))~{})~{}))"
))

df_2s_unp <- data.frame(equations = c(
  "A==N(mu[A]==1,sigma[A]==1)",
  "B==N(mu[B]==1,sigma[B]==1)",
  "d==N(mu[A]-mu[B],sigma[A]^2+sigma[B]^2)",
  "MXL==max*({}~abs({}~italic(bar(d))%+-%1.96~s[diff]~{})~{})",
  "Rel~MXL==frac(max*({}~abs({}~italic(bar(d))%+-%1.96~s[diff]~{})~{}),min*({}~abs({}~italic(bar(d))%+-%italic(t)[list(0.975,n-1)]~frac(s[diff],sqrt(n))~{})~{}))"
))

tt = ttheme_minimal(core=list(fg_params=list(hjust=0,x=0.02,fontsize=10,fontfamily="serif",parse=TRUE)))



# One sample test
a=ggplotGrob(p1_1)
gs[[6]] <- p1_1
gs[[7]] <- p1_2
gs[[8]] <- tableGrob( d = df_1s,rows=NULL, cols=NULL,theme=tt)
# Two sample, paired
gs[[10]] <- p2_1
gs[[11]] <- p2_2
gs[[12]] <- tableGrob( d = df_2s_p,rows=NULL, cols=NULL,theme=tt)
# Two sample, unpaired
gs[[14]] <- p3_1
gs[[15]] <- p3_2
gs[[16]] <- tableGrob( d = df_2s_unp,rows=NULL, cols=NULL,theme=tt)


# Arrange grob obejcts into grid
gt <- arrangeGrob(grobs = gs, layout_matrix = rbind(seq(1,4,1),
                                                    seq(5,8,1),
                                                    seq(9,12,1),
                                                    seq(13,16,1)))

# Add outlines around columna and row names
gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                              t = 1, b = 1, l = 2, r = ncol(gt))

gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                              t = 2, b = nrow(gt), l = 1, r = 1)

# Change relative height of gtable rows
gt$heights <- unit(c(.04, .32,.32, .32), "npc")
gt$widths <- unit(c(.035, .325, .325, .3), "npc")




gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))




ggsave("figure/figure_1_maxel_definition.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 7, units = "in",
       dpi = 300, limitsize = TRUE,paper="letter")




# 
# 
# 
# # Table of MAXEL equations
# equ_tbl = data.frame(maxel = c("d == bar(X[A]) + bar(X[B])",
#                                "d == bar(X[A]) + bar(X[B])",
#                                "d == bar(X[A]) + bar(X[B])"), 
#                      rel_maxel = c("d==bar(X[A])+bar(X[B])", 
#                                    "d == bar(X[A]) + bar(X[B])",
#                                    "d == bar(X[A]) + bar(X[B])"))
# colnames(equ_tbl) <- c("MAXEL", "Rel. MAXEL")
# tt <- ttheme_default(core=list(fg_params = list(parse=TRUE)))
# g <- tableGrob(equ_tbl, theme=tt)
# g$heights <- unit(rep(1/nrow(g), nrow(g)), "npc")
# grid.draw(g)