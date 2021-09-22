
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











## Explore trends between MDM and CI at a population level
#
#--------------------------------------------------------------------------------------#

# Generate 1000 samples, loop through different shifts, and quantify MDM, UCL_95, UCL_90
mus = seq(-2,2,.25)
sigma = 1
n_samples = 1000
n_obs = 6
set.seed(1)
# Sample around mean


# Make list of data frames to be concatenated at end of computations
df_list <- list()
for (n in seq(1,length(mus),1)) {
  
  y_sweep <- matrix(rnorm(n_samples*n_obs,mus[n],sigma), nrow = n_samples, byrow = TRUE)
  # Calculate MDM and max abs confidence intervals
  mdm_95    <- apply(y_sweep, 1, mdm_tdist)
  # Two tailed confidence intervals
  macb_95   <- apply(y_sweep, 1, function (x)  
    macb_tdist_2sample(x=x, y=NULL, conf.level = .95))
  macb_975  <- apply(y_sweep, 1, function (x)  
    macb_tdist_2sample(x=x, y=NULL, conf.level = .975))

  mdm_diff <- mdm_95 - macb_95
  ci_diff <- macb_975 - macb_95
  coeff_mdm95 <- mdm_diff/ci_diff
  print(mus[n])
  df_list[[n]] = tibble(n = as.factor(n), mu = mus[n],
                        coeff_mdm95 = coeff_mdm95, mdm_95 = mdm_95, 
                        macb_95 = macb_95, macb_975 = macb_975)
  df_list[[n]]$mu <- as.factor(mus[n])
}

df <- do.call(rbind, df_list)     #ldply(df_list, rbind)
df$mu <- as.factor(df$mu)

# df$mu <- as.factor(df$mu)
# Get groups means and CI of mean
df_plotted <- df %>% group_by(mu) %>% 
  summarize(mean_coeff_mdm95 = mean(coeff_mdm95)) 
# lcl_mu = boot.ci(boot(df$coeff_mdm95, statistic=
#                         function(data, i) {return(mean(data[i]))}, R=1000), 
#                  conf=0.95, type="bca")$bca[4],
# ucl_mu = boot.ci(boot(df$coeff_mdm95, statistic=
#                         function(data, i) {return(mean(data[i]))}, R=1000), 
#                  conf=0.95, type="bca")$bca[5])
df_plotted$mu <- as.factor(df_plotted$mu)
# Pairwise test between groups
pw_means <- pairwise.t.test(df$coeff_mdm95, df$mu, p.adjust.method = "bonferroni",
                            paired = TRUE, alternative = "two.sided")
adj_sig_str <- adjacent_compare_str(pw_means$p.value<0.05,'*')

# Plotting MDM vs coefficient population level
gg <- ggplot(df, aes(x = mu, y=coeff_mdm95)) +
  geom_violin(scale = "width", fill = "grey85", color="white") +
  geom_point(data=df_plotted, aes(x=mu,y=mean_coeff_mdm95), color ="#FF7F00") +
  geom_hline(aes(yintercept = 0), color = "#377eb8",  size = 1, alpha = 0.2) +
  geom_hline(aes(yintercept = 1), color = "#e41a1c", size = 1, alpha = 0.2) +
  xlab(expression(mu[DM]/sigma[DM])) + ylab(expression("Coeff."~"95%"~delta[M])) +
  # geom_text(data = data.frame(), aes(x = as.factor(mus), y = rep(1.03, length(mus)), 
  #                                    label=adj_sig_str),  size=5, position = position_nudge(x = 0.5, y = 0), color="black")+
  coord_cartesian(ylim=c(0,1.05))+
  theme_classic(base_size=8)
gg
# Export figure to disk
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "5a_coeff_MDM_from_pop.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 5, dpi = 600) # paper="letter"
graphics.off()




