
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
library(colorspace)

setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

# choose colors for plotting
color_pal = brewer.pal(3, "Set1")
color_mxl_bounds =  brewer.pal(10,"PRGn")

# Initial parameters
nsamples = 10
ngroups = 2
set.seed(0)
dist_text_size=3;
title_font_size=9;
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
  annotate("text",x=max(data1_1$x)*.95, y=c(max(data1_1$y), max(data1_1$y)-.04)-.02,
           label = c("A[phantom(0)]","C[0]"), parse=TRUE,color=color_pal[c(2,1)],
           size = dist_text_size,hjust = c(1,1)) + 
  geom_vline(xintercept=C0, color = color_pal[1],size=.5) + 
  theme(plot.title = element_text(hjust = 1,size = title_font_size),legend.position = "none") +
  ggtitle(parse(text = paste0("list(A==N(mu[A],sigma[A]),~~C[0]==1)"))) 
p1_1

### 1.2: difference with one-sample
data1_2 <- data1_1
data1_2$x <- data1_1$x-C0

spline_a_2_2 <- smooth.spline(data1_2$x,dnorm(data1_2$x, mean = mu-C0, sd = 1, log = FALSE))
ttest_a_1_2 = tidy(t.test(x=rnorm(n=nsamples, mean = mu-C0, sd = sd[1]), y=NULL,
                        alternative ="two.sided", mu = 0, paired = FALSE, 
                        var.equal = FALSE, conf.level = 0.95))
df_mxl_range_1_2 = data.frame(x1 = 0, x2 = ttest_a_1_2$conf.high, y1 = .5, y2 = .5)
p1_2 <- ggplot(data1_2, aes(x = x, y = y, group = group, fill = color_pal[3])) +
  ggtitle(parse(text = paste0("MXL==max(({}~list(abs({}~UCL[D]~{}),
                              abs({}~LCL[D]~{}))))"))) +
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif")) +
  annotate("text",x=max(data1_2$x), .95*df_mxl_range_1_2$y2, label = "D",
           parse=TRUE,color=color_pal[3],size = dist_text_size,hjust = 1) + 
  # Draw CI bounds above plots
  geom_segment(aes(x = ttest_a_1_2$conf.low, y = predict(spline_a_2_2, ttest_a_1_2$conf.low)$y, 
                   xend = ttest_a_1_2$conf.low, yend = y2, color = "segment"),
               data = df_mxl_range_1_2,inherit.aes = F, size = 1, color = lighten(color_pal[3],0.5)) +
  geom_segment(aes(x = ttest_a_1_2$conf.high, y = predict(spline_a_2_2, ttest_a_1_2$conf.high)$y, 
                   xend = ttest_a_1_2$conf.high, yend = y2, color = "segment"), 
               data = df_mxl_range_1_2,inherit.aes = F, size = 1, color = lighten(color_pal[3],0.5)) +
  # Draw bracketed range for MXL
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), 
             data = df_mxl_range_1_2,inherit.aes = F, size=1, color="black", lineend="square") +
  geom_segment(aes(x = x1, y = y1-.01, xend = x1, yend = y2, colour = "segment"), 
               data = df_mxl_range_1_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_segment(aes(x = x2, y = y1-.01, xend = x2, yend = y2, colour = "segment"), 
               data = df_mxl_range_1_2,inherit.aes = F, size=1, color="black", lineend="square") +
  geom_ribbon(data = data1_2, aes(x = x,ymax = y),ymin = 0, alpha=0.3, color = color_pal[3], 
             fill = color_pal[3])
#p1_2


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
  annotate("text",x=max(data2_1$x)*.9, y=c(max(data2_1$y), max(data1_1$y)-.04)-.02,
           label = c("A","B"), parse=TRUE,color=color_pal[c(1,2)],
           size = dist_text_size,hjust = c(.5,.5)) + 
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none") +
  ggtitle(parse(text = paste0("list(A==N(mu[A],sigma[A]),B==N(mu[B],sigma[B]))"))) 
#p2_1

### 2.2: difference two sample paired
x2_2 <- seq(-3,6,.01)
spline_a_2_2 <- smooth.spline(x2_2,dnorm(x2_2, mean = mu[1], sd = 1, log = FALSE))
data2_2 = tibble(group=factor("1"), x=x2_2, y=spline_a_2_2$y)
ttest_a_2_2 = tidy(t.test(x=rnorm(n=nsamples, mean = mu[2]-mu[1], sd = 1), y=NULL,
                        alternative ="two.sided", mu = 0, paired = FALSE, 
                        var.equal = FALSE, conf.level = 0.95))
df_mxl_range_2_2 = data.frame(x1 = 0, x2 = ttest_a_2_2$conf.high, y1 = 0.5, y2 = 0.5)
p2_2 <-ggplot(data2_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  ggtitle(parse(text = paste0("MXL==max(({}~list(abs({}~UCL[D]~{}),abs({}~LCL[D]~{}))))")))+
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif")) +
  annotate("text",x=max(data2_2$x), .95*df_mxl_range_2_2$y2, label = "D",
           parse=TRUE,color=color_pal[3],size = dist_text_size,hjust = 1) +
  geom_segment(aes(x = ttest_a_2_2$conf.low, y = predict(spline_a_2_2, ttest_a_2_2$conf.low)$y, 
                   xend = ttest_a_2_2$conf.low, yend = y2, color = "segment"),
               data = df_mxl_range_2_2,inherit.aes = F, size = 1, color = lighten(color_pal[3],0.5)) +
  geom_segment(aes(x = ttest_a_2_2$conf.high, y = predict(spline_a_2_2, ttest_a_2_2$conf.high)$y, 
                   xend = ttest_a_2_2$conf.high, yend = y2, color = "segment"), 
               data = df_mxl_range_2_2,inherit.aes = F, size = 1, color = lighten(color_pal[3],0.5)) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), 
             data = df_mxl_range_2_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_segment(aes(x = x1, y = y1-.01, xend = x1, yend = y2, colour = "segment"), 
               data = df_mxl_range_2_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_segment(aes(x = x2, y = y1-.01, xend = x2, yend = y2, colour = "segment"), 
               data = df_mxl_range_2_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data2_2$y,aes(ymin=0), alpha = 0.3, fill=color_pal[3],color=color_pal[3])
#p2_2
  
  
### Case 3: Unpaired two-sample
# 3.1 distribution A, B
x <- seq(-3,6,.01)
mu <- c(1,1.6)
sd <- c(1,1)
data3_1 = rbind(
  tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=("2"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))
p3_1 <- ggplot(data3_1, aes(x=x, y=y, group=group, fill=group)) +
  geom_line(size=1,aes(color=group)) + 
  geom_ribbon(data=data3_1, aes(x=x,ymax=y),ymin=0,alpha=0.3) +
  annotate("text",x=max(data3_1$x)*.9, y=c(max(data3_1$y), max(data1_1$y)-.04)-.02,
           label = c("A","B"), parse=TRUE,color=color_pal[c(1,2)],
           size = dist_text_size,hjust = c(.5,.5)) +  
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none")  +
  ggtitle(parse(text = paste0("list(A==N(mu[A],sigma[A]),B==N(mu[B],sigma[B]))"))) 


### 3.2: difference of Unpaired two-sample
x3_2 <- seq(-3,6,.01)
spline_a_3_2 <- smooth.spline(x3_2,dnorm(x3_2, mean = mu[1], sd = sqrt((sd[1]^2+sd[2]^2/2)), 
                                         log = FALSE))
data3_2 = tibble(group=factor("1"), x=x3_2, y=spline_a_3_2$y)
ttest_a_3_2 = tidy(t.test(x=rnorm(n=nsamples, mean = mu[2]-mu[1], sd = 1), y=NULL,
                          alternative ="two.sided", mu = 0, paired = FALSE, 
                          var.equal = FALSE, conf.level = 0.95))
df_mxl_range_3_2 = data.frame(x1 = 0, x2 = ttest_a_3_2$conf.high, y1 = 0.5, y2 = 0.5)
p3_2 <-ggplot(data3_2, aes(x=x, y=y, group=group, fill=group, color=color_pal[3])) +
  ggtitle(parse(text = paste0("MXL==max(({}~list(abs({}~UCL[D]~{}),abs({}~LCL[D]~{}))))")))+
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif")) +
  annotate("text",x=max(data3_2$x), .95*df_mxl_range_3_2$y2, label = "D",
           parse=TRUE,color=color_pal[3],size = dist_text_size,hjust = 1) +
  geom_segment(aes(x = ttest_a_3_2$conf.low, y = predict(spline_a_3_2, ttest_a_3_2$conf.low)$y, 
                   xend = ttest_a_3_2$conf.low, yend = y2, color = "segment"), size = 1,
               data = df_mxl_range_3_2,inherit.aes = F, color = lighten(color_pal[3],0.5)) +
  geom_segment(aes(x = ttest_a_3_2$conf.high, y = predict(spline_a_3_2, ttest_a_3_2$conf.high)$y, 
                   xend = ttest_a_3_2$conf.high, yend = y2, color = "segment"), 
               data = df_mxl_range_3_2,inherit.aes = F, size = 1, color = lighten(color_pal[3],0.5)) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), 
               data = df_mxl_range_3_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_segment(aes(x = x1, y = y1-.01, xend = x1, yend = y2, colour = "segment"), 
               data = df_mxl_range_3_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_segment(aes(x = x2, y = y1-.01, xend = x2, yend = y2, colour = "segment"), 
               data = df_mxl_range_3_2,inherit.aes = F, size=1, color="black",lineend="square") +
  geom_line(size=1,color=color_pal[3]) + 
  geom_ribbon(ymax = data3_2$y,aes(ymin=0), alpha = 0.3, fill=color_pal[3],color=color_pal[3])
#p3_2


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


### Equation column in figure
df_1s <- data.frame(equations = c(
  "list(italic(bar(x)[D])==bar(x)[A]-C[0],~s[D]^2==s[A]^2)",
  "CL[D]==italic(bar(x)[D])-C[0]%+-%italic(t)[alpha/2]~s/sqrt(n)*phantom(0)*phantom(0)*phantom(0)*phantom(0)*
  phantom(0)*phantom(0)*phantom(0)*phantom(0)*phantom(0)*phantom(0)"
  ))


df_2s_p <- data.frame(equations = c(
  "italic(x[list(i,D)])==italic(x[list(i,A)]) - italic(x[list(i,B)])",
  "list(italic(bar(x)[D])==frac(sum(italic(x[list(i,D)])),2), 
  italic(s[D]^2)==frac(sum((italic(x[list(i,D)])-italic(bar(x)[D])))^2,2))",
  "CL[D]==italic(bar(x)[D]) %+-% italic(t[alpha/2])~s[D]/sqrt(n)",
  "Rel~MXL==frac(MXL,min((list(abs(~UCL[D]) , abs(~LCL[D])))))"
))
# sum(x[i], i==1, n)

df_2s_unp <- data.frame(equations = c(
  "list(italic(bar(x)[D])==italic(bar(x)[A]) - italic(bar(x)[B]),phantom(0)*italic(s[D]^2)==s[A]^2 + s[B]^2)",
  "CL[D]==italic(bar(x)[D]) %+-% italic(t[alpha/2])~s[D]/sqrt(n)",
    "Rel~MXL==frac(MXL,min((list(abs(~LCL[D]) , abs(~UCL[D])))))"
))

tt = ttheme_minimal(core=list(fg_params=list(hjust=0,x=0.04,fontsize=10,
                                             fontfamily="serif",parse=TRUE)))
# modify theme in place theme=modifyList(tt, list(core=list(fg_params=list(x=-.03))))


# One sample test
gs[[6]] <-  p1_1 + labs(tag = "A", face="bold")
gs[[7]] <-  p1_2 + labs(tag = "B", face="bold")
gs[[8]] <-  tableGrob(d = df_1s, rows=NULL, cols=NULL, theme=modifyList(tt, list(core=list(fg_params=list(x=.05)))))
# Two sample, paired
gs[[10]] <- p2_1+ labs(tag = "C", face="bold")
gs[[11]] <- p2_2+ labs(tag = "D", face="bold")
gs[[12]] <- tableGrob( d = df_2s_p, rows=NULL, cols=NULL, theme=tt)
# Two sample, unpaired
gs[[14]] <- p3_1+ labs(tag = "E", face="bold")
gs[[15]] <- p3_2+ labs(tag = "F", face="bold")
gs[[16]] <- tableGrob( d = df_2s_unp, rows=NULL, cols=NULL, theme=tt)


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

#gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
#                              t = 1, b = nrow(gt), l = 3, r = 3)


# Change relative height of gtable rows and columns
gt$heights <- unit(c(.04, .32,.32, .32), "npc")
gt$widths <- unit(c(.035, .3, .355, .3), "npc")
# Apply row and column heights to gtable
gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))


# Export figure to disk
ggsave("figure/figure_1_maxel_definition.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 7, units = "in",
       dpi = 300, limitsize = TRUE,paper="letter")


