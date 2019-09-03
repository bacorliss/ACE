
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
mu <- c(1,2)
sd <- c(1,1)
data1_1 = tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y)
C0 <- 0.05;
p1_1 <- ggplot(data1_1, aes(x = x, y = y, group = group, fill = color_pal[2])) +
  geom_line(size = 1, color = color_pal[2]) + 
  geom_ribbon(data = data1_1, aes(x = x,ymax = y),ymin = 0, alpha=0.3, color = color_pal[2], fill = color_pal[2]) +
  annotate("text", x = mu[1], y = 0.4 + 0.02, label = "A",
           parse=TRUE,color=color_pal[2],size = dist_text_size) + 
  geom_vline(xintercept=C0, color = color_pal[1]) + 
  annotate("text", x = C0 - .3, y = 0.4 + 0.02, label = "C[0]",
           parse=TRUE,color=color_pal[1],size = dist_text_size) +
  ggtitle(expr(atop(paste(A==N(mu[A]==1,sigma[A]==1),"\n", C[0]==phantom(), !!C0,"   ")))) + 
  theme(plot.title = element_text(hjust = 1,size = title_font_size),legend.position = "none") 


### 1.2: difference with one-sample
data2_1 <- data1_1
data2_1$x=data1_1$x-C0
p1_2 <- ggplot(data2_1, aes(x = x, y = y, group = group, fill = color_pal[2])) +
  geom_line(size = 1, color = color_pal[2]) + 
  geom_ribbon(data = data2_1, aes(x = x,ymax = y),ymin = 0, alpha=0.3, 
              color = color_pal[2], fill = color_pal[2]) +
  annotate("text", x = mu[1], y = 0.4 + 0.02, label = "d",
           parse=TRUE,color=color_pal[2],size = dist_text_size) + 
  geom_vline(xintercept=C0, color = color_pal[1]) + 
#  annotate("text", x = C0 - .3, y = 0.4 + 0.02, label = "C[0]",
#           parse=TRUE,color=color_pal[1],size = dist_text_size) +
  ggtitle(expr(paste(d==N(mu[A]-C[0],sigma[A]),"   "))) + 
  theme(plot.title = element_text(hjust = 1, size = title_font_size)) 


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
  ggtitle(expr(paste(A==N(mu[A]==1,sigma[A]==1),"\n",B==N(mu[B]==2,sigma[B]==1),"  "))) + 
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
  ggtitle(expr(paste(d == N(mu[A]-mu[B],{sigma^2}[A]+{sigma^2}[B]),"   "))) +
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
  ggtitle(expr(paste(A==N(mu[A]==1,sigma[A]==1),"\n",B==N(mu[B]==2,sigma[B]==1),"  "))) + 
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none") 


### 3.2: difference of Unpaired two-sample
x3_2 <- seq(-3,6,.01)
data3_2 = tibble(group=factor("1"), x=x3_2,
                 y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y)
p3_2 <-ggplot(data3_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data3_2$y,aes(ymin=0), alpha = 0.3, fill=color_pal[3],color=color_pal[3]) +
  annotate("text", x=mu[2]-mu[1], y=0.4+.03, label= 
             "d == N(mu[A]-mu[B],sigma^2[A] + sigma^2[B])",
           parse=TRUE,color=color_pal[3],size=dist_text_size) +
  geom_vline(xintercept=mu[2]-mu[1]-2*sd[1]) +
  geom_vline(xintercept=mu[2]-mu[1]+2*sd[2]) +
  ggtitle(expr(paste(d == N(mu[A]-mu[B],{sigma^2}[A]+{sigma^2}[B]),"   "))) +
  theme(plot.title = element_text(hjust = 1, size = 10)) 


### Figure creation
## Assign objects to GROB table
# Column names
gs <- lapply(1:20, function(ii) grobTree(rectGrob(gp = gpar(fill = ii, alpha = 0.5)), textGrob(ii)))
gs[[2]] <- textGrob("Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[3]] <- textGrob("Effect", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[4]] <- textGrob("MAXEL", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[5]] <- textGrob("Rel. MAXEL", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))

# Row Names for figure
gs[[1]] <- textGrob("", rot = 90)
gs[[6]] <- textGrob("1-Sample", rot = 90, just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[11]] <- textGrob("2-Sample, Paired", rot = 90, just = "centre", gp=gpar(fontface = "bold", fontsize = 10))
gs[[16]] <- textGrob("2-Sample, Unpaired", rot = 90, just = "centre", gp = gpar(fontface = "bold", fontsize = 10))

### Equations for figure
maxel_1s = expression(MXL==italic(max)(phantom(.)*
                             abs(phantom(.)*italic(bar(x)[1]) - beta %+-% 
                             1.96~s[1]*phantom(.))*phantom(.)))
maxel_2s_paired = expression(MXL==italic(max)(phantom(.)*
                                    abs(phantom(.)*italic(bar(d)) %+-% 
                                    1.96~s[diff]*phantom(.))*phantom(.)))
rmaxel_2s_paired = expression(RMXL==phantom(.)*frac(max(phantom(.)*abs(phantom(.)*italic(bar(d)) %+-% 
                              1.96~s[diff]*phantom(.))*phantom(.)),
                              min(phantom(.)*abs(phantom(.)*italic(bar(d)) %+-% 
                              italic(t)[0.975, n-1]*frac(italic(s)[diff],sqrt(n))*phantom(.)))))

# One sample test
a=ggplotGrob(p1_1)
gs[[7]] <- p1_1
gs[[8]] <- p1_2
gs[[9]] <- textGrob(maxel_1s, just = "centre",gp=gpar(fontsize=equ_font_size,fontfamily="serif"))
gs[[10]] <- textGrob(expression(""), just = "centre",gp=gpar(fontsize=equ_font_size))
# Two sample, paired
gs[[12]] <- p2_1
gs[[13]] <- p2_2
gs[[14]] <- textGrob(maxel_2s_paired, just = "centre",gp=gpar(fontsize=equ_font_size,fontfamily="serif"))
gs[[15]] <- textGrob(rmaxel_2s_paired, just = "centre",gp=gpar(fontsize=equ_font_size,fontfamily="serif"))
# Two sample, unpaired
gs[[17]] <- p3_1
gs[[18]] <- p3_2
gs[[19]] <- textGrob("d == bar(X[A]) + bar(X[B])", just = "centre",gp=gpar(fontsize=equ_font_size,fontfamily="serif"))
gs[[20]] <- textGrob("d == bar(X[A]) + bar(X[B])", just = "centre",gp=gpar(fontsize=equ_font_size,fontfamily="serif"))
# Arrange grob obejcts into grid
gt <- arrangeGrob(grobs = gs, layout_matrix = rbind(seq(1,5,1),
                                                    seq(6,10,1),
                                                    seq(11,15,1),
                                                    seq(16,20,1)))

gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                              t = 1, b = 1, l = 2, r = ncol(gt))

gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                              t = 2, b = nrow(gt), l = 1, r = 1)

# Change relative height of gtable rows
gt$heights <- unit(c(.03, .32,.32, .32), "npc")
gt$widths <- unit(c(.035, .275, .275, .2, .2), "npc")

gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))


#pdf()

ggsave("figure/figure_1_maxel_definition.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 6.5, height = 7, units = "in",
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