


source("R/mhd.R")
library(broom)

x <- rnorm(10, mean = 0, sd = 2)

# Mean Shift to zero
x <- x - mean(x)


# Test base case of centrally aligned sample, where normal CI and UCL_
ucl <- t_dist_1s(x)
t_mean_ucl = tidy(t.test(x))
print(sprintf("MHD UCL: %.3f,  Ttest CI upper: %.3f", ucl, t_mean_ucl$conf.high))


# Test base case of centrally aligned sample, where normal CI and UCL_
ucl <- t_dist_1s(x+1)
t_mean_uci = tidy(t.test(x+1, conf.level = 0.95))
t_mean_uci = tidy(t.test(x+1,, conf.level = 0.975))

print(sprintf("MHD UCL: %.3f,  Ttest CI upper: %.3f", ucl, t_mean_uci$conf.high))
