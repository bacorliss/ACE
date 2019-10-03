


source("R/mhd.R")
library(broom)

# Generate repeatable sample of random normal numbers
mu = 0
std = 1
set.seed(0)
y <- rnorm(10, mean = mu, sd = std)

# Mean Shift samples to zero
y <- y - mean(y)



# Define sweep for shift in sample mean from zero
x <- seq(from = 0, to = 2*std, by = 0.1)

# Replicate samples so each of the x offsets can be applied to a seprate row
y_expanded <- matrix(y, nrow = length(x),ncol = length(y), byrow = TRUE)
y_expanded + dim(x[row(y_expanded)])
y_sweep <- sweep(y_expanded,1, x,'+')

# Calculate conifdence intervals
ucl_95    <- apply(y_sweep, 1, mhd_1sample) 
uci_95  <- apply(y_sweep, 1, function (x)  mean(x) + qt(0.950, df = length(x)-1) 
                    * sd(x)/sqrt(length(x)))
uci_975  <- apply(y_sweep, 1, function (x)  mean(x) + qt(0.975, df = length(x)-1) 
                    * sd(x)/sqrt(length(x)))


ucl_unf_975  <- apply(y_sweep, 1, function(x) ucl_mean_tdist(x_bar = mean(x), sx = sd(x), 
                                     nx = length(x), dfx = length(x)-1, alpha = 0.05))



  plot(x, ucl_95)
lines(x, uci_95) 
lines(x, uci_975)


