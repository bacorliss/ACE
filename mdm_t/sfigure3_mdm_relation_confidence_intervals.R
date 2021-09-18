
#' Characterize relation between mdm and two-tailed confidence intervals.
#' Leverage relation to produce a look-up table to calculate the MDM which 
#' accelerates computation.

# Load required packages
#-------------------------------------------------------------------------------
# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(broom)
p_load(scales)
p_load(ggplot2)
p_load(dplyr)
# p_load(plyr)
p_load(grid)
p_load(gridExtra)
p_load(colorspace)
p_load("RColorBrewer")
p_load(cowplot)
p_load(boot)
p_load(tidyr)
p_load(rbenchmark)
# User defined libraries
source("R/mdm.R")
source("R/row_stats_toolbox.R")
source("R/norm_versus_fnorm.R")


# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
fig_num = "3"
dir.create(file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep="")), showWarnings = FALSE,recursive = TRUE)
rand.seed=0

# Colors for figure
# Orange: MDM, [255, 127, 0], #FF7F00
# CL90 Blue, [255, 151, 151], #FF9797
# CL95 Red, [147, 173, 221], #93ADDD
# rTOST95, Purple, [112, 48, 160], #7030A0
col_pal <- data.frame(mdm95 = "#FF7F00", CL95 ="#FFA7A7", CL90 = "#BFCFEB", rTOST95 = "#7030A0")

ztest <- function(x, y = NULL, mu = 0, conf.level= 0.95, alternative = "two.sided") {
  #' Classic Z test of one or two samples
  #' 
  #' 
  #' 
  #' 
  if (is.null(y)) { # One sample
    n_x <- length(x)
    z_stat <- (mean(x) - mu) / (sd(x)/sqrt(n_x))
  } else {          # Two sample
    n_x <- length(x); n_y = length(y)
    z_stat <- (mean(y) - mean(x)- mu) / sqrt(sd(x)^2/n_x - sd(y)^2/n_y)
  }
  
  two_tail_p   <- 2* pnorm(abs(z_stat),lower.tail = FALSE, log.p = FALSE)
  lesser_p <- pnorm(z_stat,lower.tail = TRUE, log.p = FALSE)
  greater_p <- 1 - pnorm(z_stat,lower.tail = TRUE, log.p = FALSE)
  
  if (alternative=="two.sided") {# Two_tail (not equal to): H0 u = u0, H1: u != u0
    p <- two_tail_p
  } else if (alternative=="lesser") {# Lower Tail (less than): H0 u >= u0, H1: u < u0
    p <- lesser_p
  } else if (alternative=="greater") {# Upper Tail (greater than): H0 u <= u0, H1: u > u0
    p <- greater_p 
  } else {errorCondition("ztest: unknown alternative hypothesis")}
  
  # df <- tibble(two_tail_p=two_tail_p, lesser_p = lesser_p, greater_p = greater_p)
  return(p)
}

ztost <- function(x, delta) {
  #' Two one sided Tests for z distribution, used for equivalence testing
  #' 
  #' 
  # Equivalence: -delta <=  u1 < delta
  # z_low  <- (mean(x) - delta)/(sd(x) / sqrt(length(x)))
  # z_high <- (mean(x) + delta)/(sd(x) / sqrt(length(x)))
  
  p_lower <- ztest(x, mu = -abs(delta), alternative = "greater")
  p_upper <- ztest(x, mu =  abs(delta), alternative = "lesser" )
  # H01: delta <= -delta_0,     Ha1: delta >  -delta_0
  # p1 <- pnorm(-zp_lower_low)
  # H02: delta >= delta_0,      Ha2: delta <  delta_0
  # p2 <- pnorm(z_high)
  
  p = max(c(p_lower,p_upper))
  return(p)
}

rev_ztost <- function(xs, alpha) {
  #' Calculates the reverse two one sided tests for z distribution, returns the 
  #' largest equivalence region where results are still significant (p=0.05)
  rtost_crit <- uniroot(function(x) ztost(xs, delta=x) - alpha,
                        interval = c(0, abs(mean(xs)) + 6*sd(xs) ), 
                        tol = .Machine$double.eps)$root
  return(rtost_crit)
}

RootSpline1 <- function (x, y, y0 = 0, verbose = FALSE) {
  #' Given (x, y) data, find x where the linear interpolation crosses y = y0
  #' the default value y0 = 0 implies root finding
  #' since linear interpolation is just a linear spline interpolation
  #' the function is named RootSpline1
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



# Compare MDM, CL95, CL90 through null region
#
#------------------------------------------------------------------------------
# Generate repeatable sample of random normal numbers
mu = 0
sigma = 1
sigma_range = 0.32
n_obs = 35
set.seed(0)
y <- rnorm(n_obs, mean = mu, sd = sigma)

# Mean Shift samples to zero
y <- y - mean(y)
# Define sweep for shift in sample mean from zero
df = tibble(x = seq(from = -sigma_range*sigma, to = sigma_range*sigma, by = 0.02),grp=1)
# Replicate samples so each of the x offsets can be applied to a seprate row
y_expanded <- matrix(y, nrow = length(df$x),ncol = length(y), byrow = TRUE)
y_expanded + dim(df$x[row(y_expanded)])
y_sweep <- sweep(y_expanded,1, df$x,'+')

# Calculate most hidden difference
df$mdm_95    <- apply(y_sweep, 1, mdm_tdist)
# Max absolute one tailed confidence bounds
df$macb_95   <- apply(y_sweep, 1, function (x)  
  macb_tdist_2sample(x=x, y=NULL, conf.level = .95))
df$macb_97p5  <- apply(y_sweep, 1, function (x)  
  macb_tdist_2sample(x=x, y=NULL, conf.level = .975))


# T test p value
df$coeff_mdm_95 = (df$mdm_95  - df$macb_95) / (df$macb_97p5 - df$macb_95)
df$ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
# df$mdm_95 - df$macb_95 
x_critical <-  RootSpline1(x=df$x,y=df$ttest_p_val,y0 = 0.05)
ci_labels = c(bquote(max((~abs(CI[95])))~" "),
              bquote(max((~abs(CI[90])))), bquote(MDM[95]))

## Subplot A,B: MDM[a] transitions from CI[a] to CI[a/2] as the sample mean departs from zero 
g1A = ggplot(data=df,mapping = aes(x=x,y=mdm_95))
g1A <- g1A + theme_classic() +
  geom_rect(aes(xmin=-Inf, xmax=x_critical[1], ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_rect(aes(xmin=x_critical[2], xmax=Inf, ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_line(aes(x=x,y=macb_95, col="CL_90"), size=.8, linetype="solid", color = col_pal$CL90) +
  geom_line(aes(x=x,y=macb_97p5, col="CL_95"), size=.8, linetype="solid",color = col_pal$CL95) +
  geom_point(aes(col="MDM_95"), shape = 16,size = 1, color = col_pal$mdm95) +
  xlab("") + ylab("f(x)") +
  theme_classic(base_size=8) + theme(
    #axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
    legend.position = "right",legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical",
    legend.text=element_text(size=8),legend.key.size = unit(.5,"line"), legend.text.align=0,
    axis.title.y = element_text(size = 8)) 
  # scale_y_continuous(labels=function(x) sprintf("%.1f", x))
g1A 

# Plot p value of sample sweep
g1B_labels = c("Conf. Level  ","t-test","Crit. Region")
g1B = ggplot(data=df,mapping = aes(x=x,y=ttest_p_val))
g1B <- g1B +
  geom_rect(aes(xmin=-Inf, xmax = x_critical[1], ymin = -Inf, ymax = Inf, fill = lighten('black',0.9))) +
  geom_rect(aes(xmin = x_critical[2], xmax = Inf, ymin = -Inf, ymax = Inf, fill = lighten('black',0.9))) +
  geom_hline(aes(yintercept=0.05, col="Conf. Level"), color="red") + 
  geom_point(aes(col="t-test "), shape = 16, size = 0.5, color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression(bar(x)[DM]/s[DM])) + ylab("P-value") +
  # scale_y_continuous(name="p value", limits=c(-.1, 1.1),breaks=seq(0, 1, .5)) +
  # scale_shape_manual("Legend", labels=g1B_labels, values = c(2,2,NA)) +
  # scale_linetype_manual("Legend", labels = g1B_labels, values=c("dotted", "dotted",NA))+
  # scale_color_manual("Legend", labels = g1B_labels,values=c('red','black',"grey")) +
  scale_fill_manual("",values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) + 
  theme_classic(base_size=8) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = c(0,.5,1)) + 
  coord_cartesian(ylim=c(0,1))
g1B
# Use cowplot to align subplots
cg = plot_grid(g1A, g1B, label_size = 12, ncol=1,rel_heights = c(.45,.55) )
cg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "2AB_MDM_vs_CI95.tiff",sep=""),
          g1A, ncol = 1, nrow = 1, base_height = 1.4, base_width = 3, dpi = 600) # paper="letter"
graphics.off()


  
gg <- ggplot(data = df, aes(x=x, y=coeff_mdm_95)) + 
  ylab("Relative Distance") + xlab(expression(bar(x)[DM]/s[DM])) +
  geom_rect(aes(xmin=-Inf, xmax = x_critical[1], ymin = -Inf, ymax = Inf, fill = lighten('black',0.9))) +
  geom_rect(aes(xmin = x_critical[2], xmax = Inf, ymin = -Inf, ymax = Inf, fill = lighten('black',0.9))) +
  scale_fill_manual("",values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) + 
  geom_hline(yintercept=0, linetype="solid", color = col_pal$CL90, size=.8) + 
  geom_hline(yintercept=1, linetype="solid", color = col_pal$CL95, size=.8) + 
  geom_point(shape = 16,size = 1, color = col_pal$mdm95) + 
  theme_classic(base_size = 8) + theme(legend.position = "none") 
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "2C_MDM_vs_CI95.tiff",sep=""),
          gg, base_height = 1.2, base_width = 5, dpi = 600) 
graphics.off()

