


library(broom)
library(scales)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(colorspace)
library("RColorBrewer")
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
n_obs = 10
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
# plot
p1_1 = ggplot(data=df,mapping = aes(x=x,y=mhd_95))
p1_1 <- p1_1 + theme_classic() +
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
    legend.text=element_text(size=8))
p1_1

# Plot p value of sample sweep
p1_2 = ggplot(data=df,mapping = aes(x=x,y=ttest_p_val))
p1_2 <- p1_2 +
  geom_rect(aes(xmin=-Inf, xmax=x_critical[1], ymin=-Inf, ymax=Inf,fill = "Crit. Region")) +
  geom_rect(aes(xmin=x_critical[2], xmax=Inf, ymin=-Inf, ymax=Inf,fill = "Crit. Region")) +
  geom_point(aes(col="t-test "), shape=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_hline(aes(yintercept=0.05, col="Conf. Level")) + 
  xlab(expression(bar(x))) + ylab("p value") + 
  scale_shape_manual("", values = c(2,2)) +
  scale_linetype_manual("",values=c("dotted", "dotted"))+
  scale_color_manual("", labels = c("Conf. Level","t-test "),values=c('red','black'))+
  scale_fill_manual("",values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) + 
  theme(
    legend.position = "right", legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical",
    legend.text=element_text(size=8))
p1_2

# Assemble grid for plotting
gs <- lapply(1:2, function(ii) grobTree(rectGrob(gp = gpar(alpha = 0.5))))
gs[[1]] = p1_1 #+ labs(tag = expression(bold(A)))
gs[[2]] = p1_2 #+ labs(tag = expression(bold(B)))
  
# Arrange grob obejcts into grid
gt <- arrangeGrob(grobs = gs, layout_matrix = rbind(c(1),
                                                      c(2)))
# Change relative height of gtable rows and columns
gt$heights <- unit(c(.6, .4), "npc")
# gt$widths <- unit(c(.5,.5), "npc")
# Apply row and column heights to gtable
gt_pad <- gtable::gtable_add_padding(gt, unit(0, "inch"))

# Export figure to disk
ggsave("figure/figure_2AB_MHD_vs_CI95.pdf", gt_pad, device = "pdf", path = NULL,
       scale = 1, width = 7, height = 2.5, units = "in",
       dpi = 300, limitsize = TRUE,paper="letter")
  



## Explore trends between MHD and CI at a population level
#--------------------------------------------------------------------------------------#
# Generate 1000 samples, loop through different shifts, and quantify MHD, UCL_95, UCL_90
source("R/mhd.R")
mu = c(-100, 1, 0, 1, 100)
sigma = .01
n_samples = 1000
n_obs = 10
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
  
  y_sweep = y_samples+mu[5]
  
  
  

  mhd_1sample( y_sweep[1,] )
  max(abs( conf_interval_fcn(y_sweep[1,], 0.10) ))
  max(abs( conf_interval_fcn(y_sweep[1,], 0.05) )) 

  
  max(abs( conf_range_fcn(y_sweep[1,], 0.10) ))
  
  
  
  
  
  # Calcualte most hidden difference
  mhd_95    <- apply(y_sweep, 1, mhd_1sample)
  

  mcl_95   <- apply(y_sweep, 1, function (x)  max(abs( conf_interval_fcn(x, 0.10) )))
  mcl_975  <- apply(y_sweep, 1, function (x)  max(abs( conf_interval_fcn(x, 0.05) )))
  # ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
  
  
  mhd_diff <- mhd_95 - mcl_95
  ci_diff <- mcl_975 - mcl_95
  norm_mhd <- mhd_diff/ci_diff
  mcl_975 -mhd_95
  plot(norm_mhd)
}



# 
# ucl_95    <- apply(y_samples, 1, mhd_1sample)
# if (mu<0) {
#   uci_95  <- apply(y_samples, 1, function (x)  mean(x) - qt(0.950, df = length(x)-1) 
#                    * sd(x)/sqrt(length(x)))
#   uci_975  <- apply(y_samples, 1, function (x)  mean(x) - qt(0.975, df = length(x)-1) 
#                     * sd(x)/sqrt(length(x)))
# }else {
#   uci_95  <- apply(y_samples, 1, function (x)  mean(x) + qt(0.950, df = length(x)-1) 
#                          * sd(x)/sqrt(length(x)))
#   uci_975  <- apply(y_samples, 1, function (x)  mean(x) + qt(0.975, df = length(x)-1) 
#                           * sd(x)/sqrt(length(x)))
# }


ci_df = data.frame(mu=mu,ci_type = c(rep('MHD 95',n_samples), rep('UCI 95',n_samples), 
                                       rep('UCI 95.7',n_samples)), 
                     ci_value=c(ucl_95, uci_95, uci_975))


ci_df %>% group_by(ci_type) %>% summarise_each(list(mean=mean,sd=sd, 
                                                    std.err= function(x) sd(x)/sqrt(n_obs)))


p <- ggplot(data=ci_df, aes(x=ci_type, y=ci_value)) + 
  geom_errorbar(aes(ymin=mean(x)-sd(x),ymax=mean(x)+sd(x)))
p

p <- ggplot(data=ci_df, aes(x=ci_type, y=ci_value)) + 
  geom_boxplot()
p