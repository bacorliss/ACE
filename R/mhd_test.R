


library(broom)
library(scales)
library(ggplot2)
library(dplyr)
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
  #
  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]; y <- y[ind]
  }
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] <= 0)
  ## analytical root finding
  xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  ## make a plot?
  if (verbose) {
    plot(x, y, "l"); abline(h = y0, lty = 2)
    points(xr, rep.int(y0, length(xr)))
  }
  ## return roots
  xr
}


# Generate repeatable sample of random normal numbers
mu = 0
std = 1
n_obs = 10
set.seed(0)
y <- rnorm(n_obs, mean = mu, sd = std)

# Mean Shift samples to zero
y <- y - mean(y)

# Define sweep for shift in sample mean from zero
df = tibble(x = seq(from = -1*std, to = 1*std, by = 0.01))

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


x_critical <-  RootSpline1(x=df$x,y=df$ttest_p_val,y0 = 0.95)

# plot
p1_1 = ggplot(data=df,mapping = aes(x=x,y=mhd_95))
p1_1 <- p1_1 +
  geom_rect(aes(xmin=x_critical[1], xmax=x_critical[2], ymin=-Inf, ymax=Inf,fill = "Noncrit. Region")) +
  geom_point(shape=21, aes(col="MHD[alpha]")) + 
  geom_line(aes(x=x,y=mcl_95 , col="Max(|CI[2*alpha]|)"), linetype=2, lwd = 1) +
  geom_line(aes(x=x,y=mcl_975, col="Max(|CI[2*alpha]|)"), linetype=1, lwd = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression('bar(x)')) + ylab("f(x)") +
  scale_linetype_manual(values=c("twodash", "dotted")) +
  scale_color_manual(values=c('blue','blue','black')) +  
  scale_fill_manual("", values = lighten('red',.8), guide = guide_legend(override.aes = list(alpha = 1))) 
  #scale_color_manual(labels = c("MHD[alpha]", "Max(|CI[2*alpha]|)", "Max(|CI[alpha]|)"))
p1_1
legend(1, 1, legend=c("MHD_95", "UCL_95", "UCL_90", "- UCL_95", "- UCL_90"),
       col=c("black", alpha(c("blue","red","blue","red"),0.25)), lty=c(NA, 1, 1, 2, 2), 
       lwd=c(1,2,2,2,2),cex=0.8,pch = c(1, NA, NA, NA, NA))

# Plot p value of sample sweep
p1_2 = ggplot(data=df,mapping = aes(x=x,y=ttest_p_val))
p1_2 <- p1_2 +
  geom_rect(aes(xmin=x_critical[1], xmax=x_critical[2], ymin=-Inf, ymax=Inf), fill = lighten('red',.8)) +
  geom_point(shape=21) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression('bar(x)')) + ylab('p value')
p1_2

# x_critical = RootSpline1(x=x,y=ttest_p_val,y0 = 0.95)
# plot(x, ttest_p_val,panel.first = rect(x_critical[1], -1e6, x_critical[2], 1e6, col=lighten('red',.8), 
#                                        border=NA), xlab=expression("bar(x)"), ylab="p value")
# lines(c(x[1], tail(x,n=1)),c(0.95,0.95),col='black',lty=2)




legend(1, 1, legend=c("MHD_95", "UCL_95", "UCL_90", "- UCL_95", "- UCL_90"),
       col=c("black", alpha(c("blue","red","blue","red"),0.25)), lty=c(NA, 1, 1, 2, 2), 
       lwd=c(1,2,2,2,2),cex=0.8,pch = c(1, NA, NA, NA, NA))




# Generate 1000 samples, loop through different shifts, and quantify MHD, UCL_95, UCL_90
mu=0
sigma=1
n_samples = 1000
n_obs=10
y_samples <- matrix(rnorm(n_samples*n_obs,mu,sigma), nrow = n_samples, byrow = TRUE)

# 
ucl_95    <- apply(y_samples, 1, mhd_1sample)
if (mu<0) {
  uci_95  <- apply(y_samples, 1, function (x)  mean(x) - qt(0.950, df = length(x)-1) 
                   * sd(x)/sqrt(length(x)))
  uci_975  <- apply(y_samples, 1, function (x)  mean(x) - qt(0.975, df = length(x)-1) 
                    * sd(x)/sqrt(length(x)))
}else {
  uci_95  <- apply(y_samples, 1, function (x)  mean(x) + qt(0.950, df = length(x)-1) 
                         * sd(x)/sqrt(length(x)))
  uci_975  <- apply(y_samples, 1, function (x)  mean(x) + qt(0.975, df = length(x)-1) 
                          * sd(x)/sqrt(length(x)))
}


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