


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
source("R/mmd.R")

fig_num = "4"
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)


# Custom functions
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


# Generate repeatable sample of random normal numbers
mu = 0
sigma = 1
sigma_range = 0.3
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



# # Most confidence limit: max(abs( confident limits() ))
# conf_interval_fcn = function(x, alpha) {
#   mean(x) + c(pnorm(1-(alpha/2), df = length(x)-1) * sd(x)/sqrt(length(x)), 
#               pnorm(   alpha/2,  df = length(x)-1) * sd(x)/sqrt(length(x)))
#   }


# Calculate most hidden difference
# source("R/mmd.R")
df$mmd_95    <- apply(y_sweep, 1, mmd_normal_zdist)
# Max absolute confidence level
df$mcl_90   <- apply(y_sweep, 1, function (x)  
  max_abs_cl_mean_z_nonstandard(mean(x), sd(x)/sqrt(length(x)), alpha=0.10) ) 
df$mcl_95  <- apply(y_sweep, 1, function (x)  max(abs( 
  max_abs_cl_mean_z_nonstandard(mean(x), sd(x)/sqrt(length(x)), alpha=0.05) )))

# T test p value
df$ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
# df$mmd_95 - df$mcl_90 
x_critical <-  RootSpline1(x=df$x,y=df$ttest_p_val,y0 = 0.05)



ci_labels = c(bquote(max((~abs(CI[95])))~" "),
              bquote(max((~abs(CI[90])))), 
              bquote(MMD[95]))
                            
## Subplot A,B: MMD[a] transitions from CI[a] to CI[a/2] as the sample mean departs from zero 
g1A = ggplot(data=df,mapping = aes(x=x,y=mmd_95))
g1A <- g1A + theme_classic() +
  geom_rect(aes(xmin=-Inf, xmax=x_critical[1], ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_rect(aes(xmin=x_critical[2], xmax=Inf, ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_line(aes(x=x,y=mcl_90, col="CI_97.5"), size=.8, linetype="dotted") +
  geom_line(aes(x=x,y=mcl_95, col="CI_95"), size=.8, linetype="dotted") +
  geom_point(aes(col="MMD_95"), shape = 1,size = 1) +
  xlab("") + ylab("f(x)") +
  scale_shape_manual("", labels=ci_labels,values=c(NA,NA,1)) +   
  scale_linetype_manual("", labels=ci_labels, values=c("solid","solid","blank")) +
  scale_color_manual("",labels=ci_labels, values = c(lighten("blue",0.6),lighten("red",0.6),"black")) +
  guides(color = guide_legend(override.aes = list(
    linetype = c("solid","solid","blank"),
    shape = c(NA,NA,1), color = c(lighten("blue",0.6),lighten("red",0.6),"black"))))   + 
  theme(
    axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
    legend.position = "right",legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical",
    legend.text=element_text(size=8),legend.key.size = unit(.5,"line"), legend.text.align=0,
    axis.title.y = element_text(size = 8)) + 
  scale_y_continuous(labels=function(x) sprintf("%.1f", x))
g1A 

# Plot p value of sample sweep
g1B_labels = c("Conf. Level  ","t-test","Crit. Region")

g1B = ggplot(data=df,mapping = aes(x=x,y=ttest_p_val))
g1B <- g1B +
  geom_rect(aes(xmin=-Inf, xmax = x_critical[1], ymin = -Inf, ymax = Inf, fill = "Crit. Region")) +
  geom_rect(aes(xmin = x_critical[2], xmax = Inf, ymin = -Inf, ymax = Inf, fill = "Crit. Region")) +
  geom_hline(aes(yintercept=0.05, col="Conf. Level")) + 
  geom_point(aes(col="t-test "), shape = 16, size = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression(bar(x))) +
  scale_y_continuous(name="p value", limits=c(-.1, 1.1),breaks=seq(0, 1, .5)) +
  scale_shape_manual("Legend", labels=g1B_labels, values = c(2,2,NA)) +
  scale_linetype_manual("Legend", labels = g1B_labels, values=c("dotted", "dotted",NA))+
  scale_color_manual("Legend", labels = g1B_labels,values=c('red','black',"grey")) +
  scale_fill_manual("",values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) + 
  theme(
    legend.position = "right", legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical",
    legend.text=element_text(size=8),legend.key.size = unit(.5,"line"),
    legend.spacing.y = unit(0, "cm"),axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 10))
g1B


# Use cowplot to align subplots
cg = plot_grid(g1A, g1B, label_size = 12, ncol=1,rel_heights = c(.5,.5) )
cg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "2AB_MMD_vs_CI95.tiff",sep=""),
          cg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 6, dpi = 600) # paper="letter"
graphics.off()



## Explore trends between MMD and CI at a population level
#--------------------------------------------------------------------------------------#

# Generate 1000 samples, loop through different shifts, and quantify MMD, UCL_95, UCL_90
mu = c(-1.0, -0.4, -0.2, 0, 0.2, 0.40, 1.0, NA)
sigma = 1
n_samples = 1000
n_obs = 35
set.seed(0)
# Sample around mean
y_samples <- matrix(rnorm(n_samples*n_obs,0,sigma), nrow = n_samples, byrow = TRUE)

# Most confidence limit: max(abs( confident limits() ))
# conf_interval_fcn = function(x, alpha) mean(x) + c(qt(1-(alpha/2), df = length(x)-1) * sd(x)/sqrt(length(x)),
#                                                    qt(alpha/2, df = length(x)-1) * sd(x)/sqrt(length(x)))

# Make list of data frames to be concatenated at end of computations
df_list <- list()
for (n in seq(1,length(mu),1)) {
  shift = mu[n]
  
  if (is.na(shift)) {
    y_sweep = y_samples - rowMeans(y_samples) + runif(n_samples, min=-.00001, max=.00001)
    shift <- "Centered"
  }else {
    y_sweep = y_samples+shift
  }
  
  # Calculate MMD and max abs confidence intervals
  mmd_95    <- apply(y_sweep, 1, mmd_normal_zdist)
  mcl_90   <- apply(y_sweep, 1, function (x)  
    max_abs_cl_mean_z_nonstandard(mean(x), sd(x)/sqrt(length(x)), alpha=0.10) )
  
  mcl_95  <- apply(y_sweep, 1, function (x)  
    max_abs_cl_mean_z_nonstandard(mean(x), sd(x)/sqrt(length(x)), alpha=0.05) )
  
  # mcl_90_t   <- apply(y_sweep, 1, function (x)  max(abs( 
  #   t.test(x, conf.level = 1-0.10)$conf.int )))
  # mcl_95_t  <- apply(y_sweep, 1, function (x)  max(abs( 
  #   t.test(x, conf.level = 1-0.05)$conf.int )))
  
  
  # ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
  
  mmd_diff <- mmd_95 - mcl_90
  ci_diff <- mcl_95 - mcl_90
  fract_mmd_95 <- mmd_diff/ci_diff
  
  df_list[[n]] = tibble(n=as.factor(n), mu = as.factor(shift), 
                        fract_mmd_95 = fract_mmd_95, mmd_95 = mmd_95, 
                        mcl_90 = mcl_90, mcl_95 = mcl_95)
}

df <- ldply(df_list, rbind)

# Get groups means and CI of mean
df_plotted <- df %>% group_by(mu) %>% 
  summarize(mean_fract_mmd_95 = mean(fract_mmd_95), 
            lcl_mu = mean(fract_mmd_95) + qnorm(0.05/(2*choose(length(mu),2))) * 
              sd(fract_mmd_95) / length(fract_mmd_95),
            ucl_mu = mean(fract_mmd_95) - qnorm(0.05/(2*choose(length(mu),2))) * 
              sd(fract_mmd_95) / length(fract_mmd_95))
df_plotted$unique_sig = as.factor(rep("#",length(mu)))

ptest_result <- pairwise.t.test(df$fract_mmd_95, df$mu, p.adjust.method = "bonferroni",
                                paired = TRUE, alternative = "two.sided")



hospital_labeller <- function(variable,ind){
  return(df_plotted$unique_sig[value])
}


unique_sig = as.factor(mu)
names(unique_sig) <- as.factor(rep("#",length(mu)))

# Plotting
g1C <- ggplot(df, aes(x=mu, y=fract_mmd_95)) + 
  geom_hline(aes(yintercept=1, col="CI_90"), linetype="dotted", size=0.8) + 
  geom_hline(aes(yintercept=0, col="CI_95"), linetype="dotted", size=0.8) +
  geom_violin(fill="grey", alpha= 0.7, color = NA, lwd=0) + 
  geom_point(data = df_plotted, aes(x=mu, y=mean_fract_mmd_95))+
  geom_linerange(data=df_plotted, aes(x=mu, y=mean_fract_mmd_95, ymin=lcl_mu, ymax=ucl_mu)) +
  theme_minimal() +
  facet_grid(.~mu, scales = "free", switch = "y") + 
  theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.text.align=0,
       strip.text.x = labeller(mu=rep("#",length(mu))),
       axis.title.y = element_text(size = 8), axis.title.x = element_text(size = 10),
       legend.key.size = unit(.5,"line"), legend.spacing.y = unit(0, "cm"),
       legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm")) +
 ylab('Norm. Units') + xlab(expression(mu)) +
  scale_color_manual("", labels=c( expression(max((~abs(CI[95])))), expression(max((~abs(CI[90]))))), 
                     values=c(lighten("blue",0.4), lighten("red",0.4))) + 
  guides(color = guide_legend(override.aes = list(
    linetype = c("solid","solid"))))
#g1C
gg1C <- ggplotGrob(g1C)
gg1C$widths[18] = 6*gg1C$widths[18]
grid.draw(gg1C)

# Export figure to disk
save_plot(paste("figure/F", fig_num, "/F", fig_num, "2C_MMD_vs_CI95.tiff",sep=""),
          gg1C, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 6, dpi = 600) # paper="letter"
graphics.off()
# 
# sample_data %>%
#   group_by(Location) %>%                       
#   summarise(res = list(tidy(t.test(temp, mu=35)))) %>%
#   #unnest()

# library(tidyr)
# library(broom)
# dt_res <- df %>%
#   group_by(mu) %>%                      
#   summarise_each(funs(mean, sd, p_val_0 = t.test(fract_mmd_95,mu=0,conf.level = 1-0.05/(2*length(mu)))$p.value,
#         p_val_1 = t.test(fract_mmd_95,mu=1,conf.level = 1-0.05/(2*length(mu)))$p.value),fract_mmd_95) 




