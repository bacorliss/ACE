

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



