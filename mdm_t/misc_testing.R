

xa = rnorm(1000, mean=1,sd=1)
xb = rnorm(1000, mean=10,sd=5)
t.test(xb)
x = xa / xb

y <- cumsum(x) / seq_along(x) 

y <- cumsum(x) / seq_along(x) 


plot(y)






n_samples = 1e4
n_trials = 1e4


df = data.frame(mean_xa = rep(0,n_trials), sd_xa = rep(0,n_trials))

for (n in seq(1, n_trials,1)) {
  xa = rnorm(n_samples, mean=1,sd=3)/rnorm(n_samples, mean=1,sd=3)
  
  df$mean_xa[n] = mean(xa)
  df$sd_xa[n]   = sd(xa)
}

sum(sd(df$mean_xa) < df$sd_xa) / n_trials






source("R/mdm.R")
rand_seed <- 0
na1=15; nb1=15; na2=15; nb2=15;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)

# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*1 + 1
xb1 <- (b1-mean(b1))/sd(b1)*1 + 1
xa2 <- (a2-mean(a2))/sd(a2)*1 + 1
xb2 <- (b2-mean(b2))/sd(b2)*1 + 1





source("R/mdm.R")
mdm_tdist(xa1, xb1)
mdm_tdist_old(xa1, xb1)


