


library(broom)
library(scales)
library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(colorspace)
library("RColorBrewer")
library(cowplot)
clr_set = brewer.pal(n = 8, name = "Set1")
# Calculate conifdence intervals
source("R/mhd.R")

# Custom functions
RootSpline1 <- function (x, y, y0 = 0, verbose = FALSE) {
  ## Custom functions## given (x, y) data, find x where the linear interpolation crosses y = y0
  ## the default value y0 = 0 implies root finding
  ## since linear interpolation is just a linear spline interpolation
  ## the function is named RootSpline1
  if (is.unsorted(x)) { ind <- order(x); x <- x[ind]; y <- y[ind] }
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] <= 0)
  ## analytical root finding
  xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  ## make a plot?
  if (verbose) {plot(x, y, "l"); abline(h = y0, lty = 2); points(xr, rep.int(y0, length(xr)))
  }
  xr
}


# Generate repeatable sample of random normal numbers
mu = 0
sigma = 1
sigma_range = 1.2
n_obs = 5
set.seed(0)
y <- rnorm(n_obs, mean = mu, sd = sigma)

# Mean Shift samples to zero
y <- y - mean(y)

# Define sweep for shift in sample mean from zero
df = tibble(x = seq(from = -sigma_range*sigma, to = sigma_range*sigma, by = 0.1),grp=1)

# Replicate samples so each of the x offsets can be applied to a seprate row
y_expanded <- matrix(y, nrow = length(df$x),ncol = length(y), byrow = TRUE)
y_expanded + dim(df$x[row(y_expanded)])
y_sweep <- sweep(y_expanded,1, df$x,'+')

# Calcualte most hidden difference
df$mhd_95    <- apply(y_sweep, 1, mhd_1sample)

# Most confidence limit: max(abs( confident limits() ))
conf_interval_fcn = function(x, alpha) mean(x) + c(qt(1-(alpha/2), df = length(x)-1) * sd(x)/sqrt(length(x)), 
                                              qt(alpha/2, df = length(x)-1) * sd(x)/sqrt(length(x)))
df$mcl_95   <- apply(y_sweep, 1, function (x)  max(abs( conf_interval_fcn(x, 0.10) ))) 
df$mcl_975  <- apply(y_sweep, 1, function (x)  max(abs( conf_interval_fcn(x, 0.05) ))) 
df$ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
# df$mhd_95 - df$mcl_95 

x_critical <-  RootSpline1(x=df$x,y=df$ttest_p_val,y0 = 0.05)

## Subplot A,B: MHD[a] transitions from CI[a] to CI[a/2] as the sample mean departs from zero 
g1A = ggplot(data=df,mapping = aes(x=x,y=mhd_95))
g1A <- g1A + theme_classic() +
  geom_rect(aes(xmin=-Inf, xmax=x_critical[1], ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_rect(aes(xmin=x_critical[2], xmax=Inf, ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_line(aes(x=x,y=mcl_95, linetype = "CI_97.5"), col=lighten("blue",0.6), size=.8) +
  geom_line(aes(x=x,y=mcl_975, linetype = "CI_95"), col=lighten("blue",0.6), size=.8) +
  geom_point( aes(shape="MHD_95"), col="black",size = 1) +
  xlab("") + ylab("f(x)") +
  scale_shape_manual("", values = 1) +
  scale_linetype_manual("", labels=c( expression(max((~abs(CI[97.5])))),
    expression(max((~abs(CI[95])))), expression(MHD[95]) ),
    values=c("solid", "dotted", "solid")) +  
  #scale_fill_manual("", values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) +  
  theme(
    axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
    legend.position = "right",legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical", 
    legend.text=element_text(size=8)) + 
  scale_y_continuous(labels=function(x) sprintf("%.1f", x))
g1A

# Plot p value of sample sweep
g1B = ggplot(data=df,mapping = aes(x=x,y=ttest_p_val))
g1B <- g1B +
  geom_rect(aes(xmin=-Inf, xmax=x_critical[1], ymin=-Inf, ymax=Inf,fill = "Crit. Region")) +
  geom_rect(aes(xmin=x_critical[2], xmax=Inf, ymin=-Inf, ymax=Inf,fill = "Crit. Region")) +
  geom_hline(aes(yintercept=0.05, col="Conf. Level")) + 
  geom_point(aes(col="t-test "), shape = 16, size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression(bar(x))) +
  scale_y_continuous(name="p value", limits=c(0, 1.1),breaks=seq(0, 1, .5)) +
  scale_shape_manual("", values = c(2,2)) +
  scale_linetype_manual("",values=c("dotted", "dotted"))+
  scale_color_manual("", labels = c("Conf. Level   ","t-test "),values=c('red','black'))+
  scale_fill_manual("",values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) + 
  theme(
    legend.position = "right", legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical",
    legend.text=element_text(size=8))
g1B



cg = plot_grid(g1A, g1B, label_size = 12, ncol=1,rel_heights = c(.5,.5) )

save_plot("figure/figure_2AB_MHD_vs_CI95.pdf", cg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 6, dpi = 300) # paper="letter"


# Assemble grid for plotting
gs <- lapply(1:2, function(ii) grobTree(rectGrob(gp = gpar(alpha = 0.5))))
gs[[1]] = g1A #+ labs(tag = expression(bold(A)))
gs[[2]] = g1B #+ labs(tag = expression(bold(B)))
  
# Arrange grob obejcts into grid
gt <- arrangeGrob(grobs = gs, layout_matrix = rbind(c(1),
                                                      c(2)))
# Change relative height of gtable rows and columns
gt$heights <- unit(c(.6, .4), "npc")
# gt$widths <- unit(c(.5,.5), "npc")
# Apply row and column heights to gtable
gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))

maxWidth = grid::unit.pmax(g1B$widths[2:5], g1B$widths[2:5])
gg1A$widths[2:5] <- as.list(maxWidth)
gg1B$widths[2:5] <- as.list(maxWidth)
grid.arrange(g1A, g1B, ncol=1)

# Export figure to disk
ggsave("figure/figure_2AB_MHD_vs_CI95.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 2.5, units = "in",
       dpi = 300, limitsize = TRUE,)
  



## Explore trends between MHD and CI at a population level
#--------------------------------------------------------------------------------------#
# Generate 1000 samples, loop through different shifts, and quantify MHD, UCL_95, UCL_90
source("R/mhd.R")
mu = c(NA, -10,-1, -.5, 0, .5, 1, 10)
sigma = 1
n_samples = 100
n_obs = 5
set.seed(0)
# Sample around mean
y_samples <- matrix(rnorm(n_samples*n_obs,0,sigma), nrow = n_samples, byrow = TRUE)

# Most confidence limit: max(abs( confident limits() ))
conf_interval_fcn = function(x, alpha) mean(x) + c(qt(1-(alpha/2), df = length(x)-1) * sd(x)/sqrt(length(x)),
                                                   qt(alpha/2, df = length(x)-1) * sd(x)/sqrt(length(x)))
conf_range_fcn = function(x, alpha)  c(qt(1-(alpha/2), df = length(x)-1) * sd(x)/sqrt(length(x)),
                                                   qt(alpha/2, df = length(x)-1) * sd(x)/sqrt(length(x)))

# Make list of dataframes to be concatenated at end of computations
df_list <- list()

for (n in seq(1,length(mu),1)) {
  
  shift = mu[n]
  
  if (is.na(shift)) {
    y_sweep = y_samples - rowMeans( y_samples)
    shift <- "Centered"
  }else {
    y_sweep = y_samples+shift
  }
  
  # Calcualte most hidden difference
  mhd_95    <- apply(y_sweep, 1, mhd_1sample)
  

  mcl_95   <- apply(y_sweep, 1, function (x)  max(abs( conf_interval_fcn(x, 0.10) )))
  mcl_975  <- apply(y_sweep, 1, function (x)  max(abs( conf_interval_fcn(x, 0.05) )))
  # ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
  
  
  mhd_diff <- mhd_95 - mcl_95
  ci_diff <- mcl_975 - mcl_95
  normalized_mhd_95 <- mhd_diff/ci_diff
  
  df_list[[n]] = tibble(n=as.factor(n), mu = as.factor(shift), normalized_mhd_95 = normalized_mhd_95, 
                        mhd_95 = mhd_95, mcl_95 = mcl_95, mcl_975 = mcl_975)
  print(df_list[[n]]$mu)
  # print()
  #print(mean(normalized_mhd_95))
}

df <- ldply(df_list, rbind)




# revalue(df, c(="Centered"))

g1C <- ggplot(df, aes(x=mu, y=normalized_mhd_95)) + 
  geom_hline(yintercept=1, col=lighten("blue",0.6), linetype = "solid", size=0.8) + 
  geom_hline(yintercept=0, col=lighten("blue",0.6), linetype = "dashed", size=0.8) +
  geom_boxplot(aes(group=n)) + 
  facet_grid(. ~ n, scales = "free",switch = "y")

gg1C <- ggplotGrob(g1C)

gg1C$widths[6] = 6*gg1C$widths[6]
grid.draw(gg1C)

# Export figure to disk
ggsave("figure/figure_2C_MHD_vs_CI95.tif", gg1C, device = "tif", path = NULL,
       scale = 1, width = 5, height = 1.5, units = "in",
       dpi = 300, limitsize = TRUE,)




sample_data %>%
  group_by(Location) %>%                       
  summarise(res = list(tidy(t.test(temp, mu=35)))) %>%
  unnest()

library(tidyr)
library(broom)
dt_res <- df %>%
  group_by(mu) %>%                      
  summarise_each(funs(mean, sd, pvalue = t.test(normalized_mhd_95,mu=0)$p.value),normalized_mhd_95) 




