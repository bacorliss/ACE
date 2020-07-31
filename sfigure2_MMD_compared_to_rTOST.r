



# Load packages
library(ggplot2)
library(tibble)
library(equivalence)
library(TOSTER)
# library(RColorBrewer)
# # library(broom)
# library(gridExtra)
# library(grid)
# library(rlang)
# library(colorspace)
# library(VGAM)
# library(boot)
library(dplyr)
library(tidyr)
# library(cowplot)

source("R/row_effect_sizes.R")
source("R/mmd.R")

ztest <- function(x, y = NULL, mu = 0, conf.level= 0.95, alternative = "two.sided") {
  
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

tost_zdist <- function(x, delta) {
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

rtost_zdist <- function(xs, alpha) {
  rtost_crit <- uniroot(function(x) tost_zdist(xs, delta=x) - alpha,
          interval = c(0, abs(mean(xs)) + 6*sd(xs) ), 
          tol = .Machine$double.eps)$root
  return(rtost_crit)
}



# Test if reverse TOST equivalent to MMD
#-------------------------------------------------------------------------------


n_sims = 100
n_samples = 1E2
n_obs = 50

mus <-  seq(-1,1,by = .05)
sigmas <- rep(.1, length(mus))

df <- tibble(mu = mus, sigma = sigmas, 
             mean_mmd95 = rep(0, length(mus)),   sd_mmd95 = rep(0, length(mus)), 
                     mean_rtost95 = length(mus), sd_rtost95 = rep(0, length(mus)),
             mean_macl90 = length(mus),          sd_macl90 = rep(0, length(mus)))

for (n in seq_along(mus)) {
  print(n)
  # Spawn samples (with params by row)
  xr <- t(sapply(1:n_samples, function(x) rnorm(n_obs, mean = mus[n], sd = sigmas[n]),
               simplify = TRUE))
  # Calculate MMD across samples
  mmd95 <- apply(xr, 1, function(x) mmd_normal_zdist(x, conf.level = 0.95))
  df$mean_mmd95[n]  <- mean(mmd95)
  df$sd_mmd95[n]    <- sd(mmd95)
  # calculate reverse TOST across all samples
  rtost95 <- apply(xr, 1, function(x) rtost_zdist(x, alpha = 0.05)  )
  df$mean_rtost95[n] <- mean(rtost95)
  df$sd_rtost95[n]   <- sd(rtost95)
  # Calculate max( |CL_95| )
  macl90 <- apply(xr, 1, function(x) max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), a = 0.10))
  df$mean_macl90[n]  <- mean(macl90)
  df$sd_macl90[n]    <- sd(macl90)
  
    
}

df_plot <- df %>% gather(metric, mean, starts_with("mean")) %>%
  gather(metric, sd, starts_with("sd"))
  


gg <- ggplot(data = df_plot,(aes(x=mu, y=mean, color = metric))) +
  geom_line() + 
  # geom_ribbon(aes(ymin=df_plot$mean-2.5*df_plot$sd, 
  #                 ymax=df_plot$mean+2.5*df_plot$sd), linetype=2, alpha=0.1) +
  theme_classic(base_size = 8)

gg
















x = rnorm(50, mean=20, sd=1); tost_p = rep(0,length(x))
delta = seq(-5,5,by =.1)
for (n in seq_along(delta)){
  tost_p[n] <- tost_zdist(x,delta[n])
  
}
plot(delta,tost_p)


# Generate samples for mu and Sd, calculate all three curves for each, the calculate SSE
mus = seq(-5,5,by=0.1)
n_samples <- length(mus)
sigmas = runif(n_samples,1,1)
n_obs = 50
salpha = 0.05

xs = t(mapply(function(x,y) rnorm(n_obs, mean = x, sd = 1), mus, sigmas))
xc = xs  - rowMeans(xs) 

inv_tost_es <- rep(0, n_samples)
mmd_es <- rep(0, n_samples)
mabs_cl <- rep(0, n_samples)
for (n in seq(1, n_samples,1)){
  
  inv_tost_es[n] <- uniroot(function(x) tost_zdist(xc[n,], delta=x) - salpha,
                            interval = c(0, abs(mean(xc[n,])) + 6*sd(xc[n,]) ), 
                            tol = .Machine$double.eps)$root
  mmd_es[n] <- mmd_normal_zdist(xc[n,])
  mabs_cl[n] <- max_abs_cl_mean_z(x_bar = mean(xc[n,]), sem_x = sd(xc[n,])/sqrt(n_obs),alpha = 0.10)
}

plot(rowSds(xc), (mmd_es - inv_tost_es))
# plot(rowSds(xc), (mmd_es-0.0445*rowSds(xc)) - inv_tost_es)



x=seq(0,10,.1)
y_sweep <- sapply(x, function(x) tost(xs[n,], y = NULL, delta = x, paired = FALSE, var.equal = FALSE,
                                      conf.level = 0.95, alpha = salpha)$tost.p.value - salpha, simplify = TRUE)
y_sweep <- sapply(x, function(x) tost_zdist(xs[1,], delta = x) - salpha, simplify = TRUE)
plot(x,y_sweep)
# 
# 
# 
# x=rnorm(50,mean=.2,sd=1)
# 
# library(BSDA)


