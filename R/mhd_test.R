


source("R/mhd.R")
library(broom)
library(scales)
library(ggplot2)
library(dplyr)

# Generate repeatable sample of random normal numbers
mu = 0
std = 1
set.seed(0)
y <- rnorm(10, mean = mu, sd = std)

# Mean Shift samples to zero
y <- y - mean(y)



# Define sweep for shift in sample mean from zero
x <- seq(from = -1.5*std, to = 1.5*std, by = 0.1)


# Replicate samples so each of the x offsets can be applied to a seprate row
y_expanded <- matrix(y, nrow = length(x),ncol = length(y), byrow = TRUE)
y_expanded + dim(x[row(y_expanded)])
y_sweep <- sweep(y_expanded,1, x,'+')

# Calculate conifdence intervals
source("R/mhd.R")
ucl_95    <- apply(y_sweep, 1, mhd_1sample)

# Right tail confidence levels
uci_right_95  <- apply(y_sweep, 1, function (x)  mean(x) + qt(0.950, df = length(x)-1) 
                    * sd(x)/sqrt(length(x)))
uci_right_975  <- apply(y_sweep, 1, function (x)  mean(x) + qt(0.975, df = length(x)-1) 
                    * sd(x)/sqrt(length(x)))
# Left tail confidence levels
uci_left_95  <- -apply(y_sweep, 1, function (x)  mean(x) - qt(0.950, df = length(x)-1) 
                       * sd(x)/sqrt(length(x)))
uci_left_975  <- -apply(y_sweep, 1, function (x)  mean(x) - qt(0.975, df = length(x)-1) 
                        * sd(x)/sqrt(length(x)))



ucl_unf_975  <- apply(y_sweep, 1, function(x) ucl_mean_tdist(x_bar = mean(x), sx = sd(x), 
                                     nx = length(x), dfx = length(x)-1, alpha = 0.05))


# Plot comparing confidence intervals to MHD
  plot(x, ucl_95, ylim = c(0, max(ucl_95)),cex=1,
       xlab=bquote("Shifted " * bar(x)), ylab='f(x)')
lines(x, uci_right_95,type="l",lty = 1,col=alpha('red', 0.25),lwd=2)
lines(x, uci_left_95,type="l",lty = 2,col=alpha('red', 0.25),lwd=2) 

lines(x, uci_right_975,lty = 1,col=alpha('blue', 0.25),lwd=2)
lines(x, uci_left_975,lty = 2,col=alpha('blue', 0.25),lwd=2)

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