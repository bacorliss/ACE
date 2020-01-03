

# This script produces a figure defining MAXEL for two-tailed upaired data in a distribution
# in a distribution free manner
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
  annotate("text",x=max(data3_1$x)*.9, y=c(max(data3_1$y), max(data1_1$y)-.03),
           label = c("A","B"), parse=TRUE,color=color_pal[c(1,2)],
           size = dist_text_size,hjust = c(.5,.5)) +  
  theme(plot.title = element_text(hjust = 1, size = title_font_size),legend.position = "none")  + 
  ylab(expression(italic(f))) +
  ggtitle(parse(text = paste0("list(A==N(mu[A],sigma[A]),B==N(mu[B],sigma[B]))")))


### 3.2: difference of Unpaired two-sample
x3_2 <- seq(0,3,.01)
spline_a_3_2 <- smooth.spline(x3_2,dnorm(x3_2, mean = mu[1], sd = 1, log = FALSE))
spline_b_3_2 <- smooth.spline(x3_2,dnorm(x3_2, mean = mu[2], sd = 1, log = FALSE))
data3_2 = rbind(
  tibble(group=factor("1"), x=x3_2, y=spline_a_3_2$y),
  tibble(group=factor("2"), x=x3_2, y=spline_b_3_2$y))
# Get confidence invterval for each distribution
ttest_a_3_2 = tidy(t.test(x=rnorm(n=nsamples, mean = mu[1], sd = 1),y=NULL, 
                          alternative ="two.sided", mu = 0, paired = FALSE, 
                          var.equal = FALSE, conf.level = 0.95))
ttest_b_3_2 = tidy(t.test(x=rnorm(n=nsamples, mean = mu[2], sd = 1),y=NULL, 
                          alternative ="two.sided", mu = 0, paired = FALSE, 
                          var.equal = FALSE, conf.level = 0.95))
# Line segment range for MAXEL
df_mxl_range_3_1 = data.frame(x1 = ttest_a_3_2$conf.low, x2 = ttest_b_3_2$conf.high, y1 = .5, y2 = .5)
# Plot
p3_2 <- ggplot(data3_2, aes(x=x, y=y, group=group, fill=group)) +
  annotate("text",x = max(data3_2$x)*.9, y = c(max(data3_2$y), max(data1_1$y)-.03), label = c("A","B"), 
           parse = TRUE,color = color_pal[c(1,2)], size = dist_text_size,hjust = c(.5,.5)) +  
  ggtitle(parse(text = paste0("MXL==max((list({}~abs({}~UCL[A]-LCL[B]~{}),
                              abs({}~LCL[A]-UCL[B]~{}))))"))) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, 
                                                            size = title_font_size,family="serif")) +
  geom_segment(aes(x = ttest_a_3_2$conf.low, y = predict(spline_a_3_2, ttest_a_3_2$conf.low)$y, 
                   xend = ttest_a_3_2$conf.low, yend = y2, color = "segment"),
               data = df_mxl_range_3_1,inherit.aes = F, size = 1, color = lighten(color_pal[1],.5)) +
  geom_segment(aes(x = ttest_a_3_2$conf.high, y = predict(spline_a_3_2, ttest_a_3_2$conf.high)$y, 
                   xend = ttest_a_3_2$conf.high, yend = y2, color = "segment"), 
               data = df_mxl_range_3_1,inherit.aes = F, size = 1, color = lighten(color_pal[1],.5)) +
  geom_segment(aes(x = ttest_b_3_2$conf.low, y =  predict(spline_b_3_2, ttest_b_3_2$conf.low)$y, 
                   xend = ttest_b_3_2$conf.low, yend = y2, color = "segment"), 
               data = df_mxl_range_3_1,inherit.aes = F, size = 1, color = lighten(color_pal[2],.5)) +
  geom_segment(aes(x = ttest_b_3_2$conf.high, y =  predict(spline_b_3_2, ttest_b_3_2$conf.high)$y, 
                   xend = ttest_b_3_2$conf.high, yend = y2, color = "segment"), 
               data = df_mxl_range_3_1,inherit.aes = F, size = 1, color = lighten(color_pal[2],.5)) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"),
               data = df_mxl_range_3_1,inherit.aes = F, size=1, color="black", lineend="square") +
  geom_segment(aes(x = x1, y = y1 - 0.02, xend = x1, yend = y2, colour = "segment"), 
               data = df_mxl_range_3_1,inherit.aes = F, size=1, color="black",linetype=1, lineend="square") +
  geom_segment(aes(x = x2, y = y1 - 0.02, xend = x2, yend = y2, colour = "segment"), 
               data = df_mxl_range_3_1,inherit.aes = F, size=1, color="black",linetype=1, lineend="square") +
  geom_line(size=1,aes(color=group)) + 
  geom_ribbon(data=data3_2, aes(x=x,ymax=y),ymin=0,alpha=0.3)  + 
  ylab(expression(italic(f)))
# p3_2

df_2s_unp <- data.frame(equations = c(
  "CL==italic(bar(x)) %+-% italic(t[alpha/2])~s/sqrt(n)",
  "Rel~MXL==frac(MXL,min((list(abs(LCL[A]) , abs(UCL[A])))))"
))


### Figure creation
## Assign objects to GROB table
# Column names
gs <- lapply(1:16, function(ii) grobTree(textGrob("")))
gs[[2]] <- textGrob("Initial Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[3]] <- textGrob("Effect Distribution", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))
gs[[4]] <- textGrob("Equations", just = "centre", gp = gpar(fontface = "bold", fontsize = 10))



# Row headers 
gs[[5]] <- grobTree(textGrob("2-S Paired", rot = 90, just = "centre", 
                             gp = gpar(fontface = "bold", fontsize = 10)),rectGrob(gp = gpar(alpha = 0.5)))


# One sample test
gs[[6]] <- p3_1 + labs(tag = "A", face="bold")
gs[[7]] <- p3_2 + labs(tag = "B", face="bold")
gs[[8]] <- tableGrob( d = df_2s_unp,rows=NULL, cols=NULL,theme=tt)



# Arrange grob obejcts into grid
gt <- arrangeGrob(grobs = gs, layout_matrix = rbind(seq(1,4,1),seq(5,8,1),seq(9,12,1),seq(13,16,1)))


# Add outlines around columna and row names
gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                              t = 1, b = 1, l = 2, r = ncol(gt))

gt <- gtable::gtable_add_grob(gt, grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                              t = 2, b = 2, l = 1, r = 1)

# Change relative height of gtable rows and columns
gt$heights <- unit(c(.04, .32,.32, .32), "npc")
gt$widths <- unit(c(.035, .3, .40, .265), "npc")
# Apply row and column heights to gtable
gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))


# Export figure to disk
ggsave("figure/figure_2_maxel_definition.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 7, units = "in",
       dpi = 300, limitsize = TRUE,paper="letter")



