
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
library(VGAM)
library(boot)

setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

# choose colors for plotting
color_pal = brewer.pal(4, "Set1")

# Initial parameters
nsamples = 50
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
sample_1_1 = rnorm(nsamples,mean=mu[1],sd=sd[1])
C0 <- 1;
p1_1 <- ggplot(data1_1, aes(x = x, y = y, group = group)) +
  geom_line(aes(color=group), size = .75) + 
  geom_ribbon(aes(x = x,ymax = y, fill=group),ymin = 0, alpha = 0.3) +
  scale_color_manual(values = color_pal[2]) + scale_fill_manual(values = color_pal[2]) +
  annotate("text",x = max(data1_1$x) * .95, y = c(max(data1_1$y), max(data1_1$y) - .04) - .02, 
           label = c("G[phantom(0)]","C[0]"), parse = TRUE,color = color_pal[c(1,2)],
           size = dist_text_size, hjust = c(1,1)) + 
  geom_vline(xintercept=C0, color = color_pal[1],size=  .75) + 
  theme(plot.title = element_text(hjust = 1,size = title_font_size),legend.position = "none") +
  ggtitle(parse(text = paste0("list(G~'~'~N(mu[G],sigma[G]),~~C[0]==1)"))) +
  ylab(expression(italic(f(x))))
p1_1

### 1.2: difference with one-sample
# Difference distribution describing effect size
ant_y1 = 0.6; 
n_dist_1_2 = data.frame(x = data1_1$x - C0, y = data1_1$y, group=factor(1))
# Folded normal distribution describing effect size
fn_dist_1_2 <- n_dist_1_2
fn_dist_1_2$y <- dfoldnorm(n_dist_1_2$x, mean = mu[1]-C0, sd = 1, log = FALSE)
# Define upper limit of confidence interveral of the mean
sample_1_2 = abs(sample_1_1-C0)#rfoldnorm(50,mean=mu-C0,sd=1, a1=1,a2=1)
ucl_1_2 <- boot.ci(boot(sample_1_2, function(d, i){mean(d[i])}, R = 10000), type = 'bca', conf=0.90)$bca[5]

p1_2 <- ggplot(n_dist_1_2, aes(x = x, y = y, group = group, fill = color_pal[3])) +
  ggtitle(parse(text = paste0("F~'~'~abs(~N(~mu[D],sigma[D])*phantom(.))"))) +
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif"),legend.position = "none") +
  annotate("text",x = max(n_dist_1_2$x) * .95, y = c(ant_y1, ant_y1 - .06) - .02, 
           label = c("F","D"), parse = TRUE, color = color_pal[c(4,3)],
           size = dist_text_size, hjust = c(1,1)) +
  geom_line(data =  n_dist_1_2, aes(x = x,y = y), size = .75, color = color_pal[3]) +
  geom_line(data = fn_dist_1_2, aes(x = x,y = y), size = .75, color = color_pal[4]) + 
  geom_ribbon(data = subset(fn_dist_1_2,x<ucl_1_2), aes(x = x,ymax = y),ymin = 0, color = color_pal[4],
              fill = color_pal[4], alpha=0.3) +
  # Draw bracketed range for MXL
  geom_segment(aes(x = 0, y = ant_y1 - 0.01, xend = 0, yend = ant_y1, colour = "segment"),
               inherit.aes = F, size=1, color = "black", lineend="square") +
  geom_segment(aes(x = 0, y = ant_y1, xend = ucl_1_2, yend = ant_y1, colour = "segment"),
             inherit.aes = F, size=1, color="black", lineend="square") +
  geom_segment(aes(x = ucl_1_2, y = ant_y1 - 0.01, xend = ucl_1_2, yend = ant_y1, colour = "segment"),
               inherit.aes = F, size=1, color="black",lineend="square") +
  ylab(expression(italic(f(x))))
p1_2



### Case 2: two sample paired
# 2.1 distribution A, B
x <- seq(-3,6,.01)
mu <- c(1,2)
sd <- c(1,1)
data2_1 = rbind(tibble(group=factor("1"), x = x, y = smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=factor("2"), x=x, y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))
sample_2_1 = rnorm(nsamples, mean=mu[2] - mu[1], sd = sd[1])
p2_1 <- ggplot(data2_1, aes(x = x, y = y, group = group)) +
  geom_line(size=1,aes(color = group)) + 
  geom_ribbon(aes(x = x,ymax = y, fill=group),ymin = 0,alpha = 0.3) +
  scale_color_manual(values = color_pal[c(1,2)]) + scale_fill_manual(values = color_pal[c(1,2)]) +
  annotate("text",x = max(data2_1$x) * .9, y = c(max(data2_1$y), max(data1_1$y) - .04) - .02,
           label = c("G","H"), parse = TRUE, color = color_pal[c(1,2)],
           size = dist_text_size, hjust = c(.5,.5)) + 
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none") +
  ggtitle(parse(text = paste0("list(G~'~'~N(mu[G],sigma[G]),H~'~'~N(mu[H],sigma[H]))"))) +
  ylab(expression(italic(f(x)))) 
p2_1

### 2.2: difference two sample paired
# Difference distribution describing effect size
ant_y2 = 0.6;
n_dist_2_2 <- data.frame(x=x,y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = 1, log = FALSE))$y,
                         group=factor(1))
# Folded normal distribution describing effect size
fn_dist_2_2 <- n_dist_2_2
fn_dist_2_2$y <- dfoldnorm(n_dist_2_2$x, mean = mu[2]-mu[1], sd = 1, a1 = 1,a2 = 1, log = FALSE)
# Simulated sample
sample_2_2 = abs(sample_2_1)
# Define upper limit of confidence interveral of the mean
ucl_2_2 <- boot.ci(boot(sample_2_2, function(d, i){mean(d[i])}, R = 10000), type = 'bca', conf=0.90)$bca[5]
# Plot difference and folded difference with UCL
p2_2 <- ggplot(n_dist_2_2, aes(x = x, y = y, fill = color_pal[3])) +
  ggtitle(parse(text = paste0("F~'~'~abs(~N(~mu[D],sigma[D])*phantom(.))")))  +
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif"),legend.position = "none") +
  annotate("text",x = max(n_dist_2_2$x) * .95, y = c(ant_y2, ant_y2 - .06) - .02, 
           label = c("F","D"), parse = TRUE,color = color_pal[c(3,4)],
           size = dist_text_size, hjust = c(1,1)) +
  geom_segment(aes(x = ucl_2_2, y = 0, xend = ucl_2_2, yend = predict(
    smooth.spline(x = fn_dist_2_2$x, y = fn_dist_2_2$y),ucl_2_2)$y, colour = "segment"),
    inherit.aes = F, size = .75, color = color_pal[3], lineend = "square") +
  geom_line(data =  n_dist_2_2, aes(x = x,y = y), size = .75,  color = color_pal[3]) +
  geom_line(data = fn_dist_2_2, aes(x = x,y = y), size = .75,color = color_pal[4]) + 
  geom_ribbon(data = subset(fn_dist_2_2,x<ucl_2_2), aes(x = x,ymax = y),ymin = 0, color = color_pal[4],
              fill = color_pal[4], alpha=0.3) +
  # Draw bracketed range for MXL
  geom_segment(aes(x = 0, y = ant_y2 - 0.01, xend = 0, yend = ant_y2, colour = "segment"),
               inherit.aes = F, size=1, color = "black", lineend="square") +
  geom_segment(aes(x = 0, y = ant_y2, xend = ucl_2_2, yend = ant_y2, colour = "segment"),
               inherit.aes = F, size=1, color="black", lineend="square") +
  geom_segment(aes(x = ucl_2_2, y = ant_y2 - 0.01, xend = ucl_2_2, yend = ant_y2, colour = "segment"),
               inherit.aes = F, size=1, color="black",lineend="square") +
  ylab(expression(italic(f(x))))
# p2_2


  
### Case 3: Unpaired two-sample
# 3.1 distribution A, B
ant_y3=0.7
x <- seq(-3,6,.01)
mu <- c(1,1.6)
sd <- c(1,1)
data3_1 = rbind(
  tibble(group=factor("1"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[1], sd = 1, log = FALSE))$y),
  tibble(group=("2"), x=x,
         y=smooth.spline(x,dnorm(x, mean = mu[2], sd = 1, log = FALSE))$y))
sample_3_1 = rnorm(nsamples, mean=mu[2] - mu[1], sd = sqrt((sd[1]^2+sd[2]^2/2)))
# Plot original distributions
p3_1 <- ggplot(data3_1, aes(x = x, y = y, group = group, fill = group)) +
  geom_line(size = 0.75,aes(color=group)) + 
  geom_ribbon(aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
  scale_color_manual(values = color_pal[c(1,2)]) + scale_fill_manual(values = color_pal[c(1,2)]) +
  annotate("text",x=max(data3_1$x)*.9, y=c(max(data3_1$y), max(data1_1$y)-.04)-.02,
           label = c("G","H"), parse=TRUE,color=color_pal[c(1,2)],
           size = dist_text_size,hjust = c(.5,.5)) +  
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none")  +
  ggtitle(parse(text = paste0("list(G~'~'~N(mu[G],sigma[G]),H~'~'~N(mu[H],sigma[H]))"))) +
  ylab(expression(italic(f(x))))
p3_1 

### 3.2: difference of Unpaired two-sample
n_dist_3_2 <- data.frame(x=x,y=smooth.spline(x,dnorm(x, mean = mu[2]-mu[1], sd = sqrt((sd[1]^2+sd[2]^2/2)),
                                                     log = FALSE))$y, group=factor(1))
# Folded normal distribution describing effect size
fn_dist_3_2 <- n_dist_3_2
fn_dist_3_2$y <- dfoldnorm(n_dist_3_2$x, mean = mu[2]-mu[1], sd = 1, a1 = 1,a2 = 1, log = FALSE)
# Simulated sample
sample_3_2 = abs(sample_3_1)#rfoldnorm(50,mean=mu-C0,sd=1, a1=1,a2=1)
# Define upper limit of confidence interveral of the mean
ucl_3_2 <- boot.ci(boot(sample_3_2, function(d, i){mean(d[i])}, R = 10000), type = 'bca', conf=0.90)$bca[5]
# Plot difference and folded difference with UCL %>% 
p3_2 <- ggplot(n_dist_3_2, aes(x = x, y = y, fill = color_pal[3])) +
  ggtitle(parse(text = paste0("F~'~'~abs(~N(~mu[D],sigma[D])*phantom(.))")))  +
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif"),legend.position = "none") +
  annotate("text",x = max(n_dist_3_2$x) * .95, y = c(ant_y3, ant_y3 - .07) - .02, 
           label = c("F","D"), parse = TRUE,color = color_pal[c(3,4)],
           size = dist_text_size, hjust = c(1,1)) +
  geom_segment(aes(x = ucl_3_2, y = 0, xend = ucl_3_2, yend = predict(
    smooth.spline(x = fn_dist_3_2$x, y = fn_dist_3_2$y),ucl_3_2)$y, colour = "segment"),
    inherit.aes = F, size = .75, color = color_pal[3], lineend = "square") +
  geom_line(data =  n_dist_3_2, aes(x = x,y = y), size = .75, color = color_pal[3]) +
  geom_line(data = fn_dist_3_2, aes(x = x,y = y), size = .75,color = color_pal[4]) + 
  geom_ribbon(data = subset(fn_dist_3_2,x<ucl_3_2), aes(x = x,ymax = y),ymin = 0, color = color_pal[4],
              fill = color_pal[4], alpha=0.3) +
  # Draw bracketed range for MXL
  geom_segment(aes(x = 0, y = ant_y3 - 0.01, xend = 0, yend = ant_y3, colour = "segment"),
               inherit.aes = F, size=1, color = "black", lineend="square") +
  geom_segment(aes(x = 0, y = ant_y3, xend = ucl_3_2, yend = ant_y3, colour = "segment"),
               inherit.aes = F, size=1, color="black", lineend="square") +
  geom_segment(aes(x = ucl_3_2, y = ant_y3 - 0.01, xend = ucl_3_2, yend = ant_y3, colour = "segment"),
               inherit.aes = F, size=1, color="black",lineend="square") +
  ylab(expression(italic(f(x))))
# p3_2


### Figure creation
## Assign objects to GROB table
# Column names
gs <- lapply(1:16, function(ii) grobTree(rectGrob(gp = gpar(fill = ii, alpha = 0.5)), textGrob(ii)))
gs[[2]] <- textGrob("Initial Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 9))
gs[[3]] <- textGrob("Effect Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 9))
gs[[4]] <- textGrob("Equations", just = "centre", gp = gpar(fontface = "bold", fontsize = 9))


# Row Names for figure
gs[[1]] <- textGrob("", rot = 90)
gs[[5]] <- textGrob("1-Sample", rot = 90, just = "centre", gp = gpar(fontface = "bold", fontsize = 9))
gs[[9]] <- textGrob("2-Sample, Paired", rot = 90, just = "centre", gp=gpar(fontface = "bold", fontsize = 9))
gs[[13]] <- textGrob("2-Sample, Unpaired", rot = 90, just = "centre", gp = gpar(fontface = "bold", fontsize = 9))


### Equation column in figure
df_1s <- data.frame(equations = c(
  "list(italic(x[list(i,D)])==italic(x[list(i,G)]-C[0]),~s[D]^2==s[G]^2)",
  'x[list(i,F)]==abs(~x[~list(i,D)]*phantom(.))',
  "UCL(F) == italic(boot)(x[list(i,F)], mean(x), 1-alpha)",
  paste0("MXL==UCL(~F~{})")
  ))


df_2s_p <- data.frame(equations = c(
  "list(italic(x[list(i,D)])==italic(x[list(i,G)]) - italic(x[list(i,H)]), ~x[list(i,F)]==abs(~x[~list(i,D)]*phantom(.)))",
  "list(italic(bar(x)[D])==frac(sum(italic(x[list(i,D)])),2), 
  ~italic(s[D]^2)==frac(sum((italic(x[list(i,D)])-italic(bar(x)[D])))^2,2))",
  "UCL(F) == italic(boot)(x[list(i,F)], mean(x), 1-alpha)",
  paste0("MXL==UCL(~F~{})")
))
# sum(x[i], i==1, n)

df_2s_unp <- data.frame(equations = c(
  "list(italic(bar(x)[D])==italic(bar(x)[G]) - italic(bar(x)[H]),phantom(0)*italic(s[D]^2)==s[G]^2 + s[H]^2)",
  "x[G]^{minute} == x[list(i,G)]-bar(x)[G]+bar(x)[GH]",
  "x[H]^{minute} == x[list(i,H)]-bar(x)[H]+bar(x)[GH]",
  "UCL(F) == italic(boot)(list(x[G]^{minute},x[H]^{minute}),abs(~diff(x)), 1-alpha)",
  paste0("MXL==UCL(~F~{})")
))

tt = ttheme_minimal(core=list(fg_params=list(hjust=0,x=0.04,fontsize=10,
                                             fontfamily="serif",parse=TRUE)))
# modify theme in place theme=modifyList(tt, list(core=list(fg_params=list(x=-.03))))


# One sample test
gs[[6]] <-  p1_1 + labs(tag = expression(bold(A)))
gs[[7]] <-  p1_2 + labs(tag = expression(bold(B)))
gs[[8]] <-  tableGrob(d = df_1s, rows=NULL, cols=NULL, 
                      theme=modifyList(tt, list(core=list(fg_params=list(x=.05)))))
# Two sample, paired
gs[[10]] <- p2_1+ labs(tag = expression(bold(C)))
gs[[11]] <- p2_2+ labs(tag = expression(bold(D)))
gs[[12]] <- tableGrob( d = df_2s_p, rows=NULL, cols=NULL, theme=tt)
# Two sample, unpaired
gs[[14]] <- p3_1+ labs(tag = expression(bold(E)))
gs[[15]] <- p3_2+ labs(tag = expression(bold(F)))
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
gt$widths <- unit(c(.035, .3, .3, .355), "npc")
# Apply row and column heights to gtable
gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))


# Export figure to disk
ggsave("figure/figure_1_maxel_definition.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 7, units = "in",
       dpi = 300, limitsize = TRUE,paper="letter")


