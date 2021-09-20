
# MDM Package: most difference in means
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.


# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(docstring)
p_load(cubature)
source("R/rationormal_toolbox.R")


# Note: when pooled variance unequal variance
# http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/poolvar.html



dft <- function(x, df, mu, sigma) {
  #' @description: Probability density function for noncentral folded t-distribution
  #' Equation from:
  #'  https://en.wikipedia.org/wiki/Folded-t_and_half-t_distributions
  #' @param x 
  #' @param df degree of freedom
  #' @param mu mean
  #' @param sigma standard deviation

  # save(list = ls(all.names = TRUE), file = "temp/dft.RData",envir = environment())
  # load(file = "temp/dft.RData")
  
  v = df
  t1 <- (1 +  (1/v) * ((x - mu)^2 / sigma^2 )) ^ ( -(v+1) / 2)
  t2 <- (1 +  (1/v) * ((x + mu)^2 / sigma^2 )) ^ ( -(v+1) / 2)
  
  d <- gamma( (v+1)/2 ) / (gamma(v/2) * sqrt(v*pi*sigma^2) ) * (t1 + t2)
  # d is zero for all x<0 by definition, so force d to 0 when x is negative
  d[x < 0] <- 0

  return(d)
}

pft <- function(upper, df, mu, sigma, lower = -Inf) {
  #' @description: cumulative distribution function for noncentral folded 
  #' t-distribution
  #' @param upper upper bounds of integration
  #' @param df degree of freedom
  #' @param mu mean
  #' @param sigma standard deviation
  #' @param lower lower bounds of integration, default -InF
  
  # save(list = ls(all.names = TRUE), file = "temp/pft.RData",envir = environment())
  # load(file = "temp/pft.RData")
  
  d <- integrate(f = function(y) dft(y, df, abs(mu),sigma),
                 lower = lower, upper = upper,
                  rel.tol = 1e-8, abs.tol = 1e-9)$value
  
  return(d)
}


qft <- function(p, df, mu, sigma, lo.bound, hi.bound) {
  #' @description Quantile function for folded t-distribution
  #'
  #' @param p probability
  # Estimate bounds to search for root

  # save(list = ls(all.names = TRUE), file = "temp/qft.RData",envir = environment())
  # load(file = "temp/qft.RData")

  # Integration is also used to calculate p, but it can be a very sparse 
  # integration (long spans where f(x)=0 and then a very small width spike in the pdf)
  # To deal with this, integration is broken up into separate intervals to make 
  # sure the spike is not missed
  
  # Value of p is always greater than mu, so integrate up to that first, but this
  # integration is split between the sparse region and region where f(x)!=0 
  # when it's near mu
  #   So we integrate 
  #   1) from 0 to max(c(lo.bound - 10*sigma,0))
  #   2) from max(c(lo.bound - 10*sigma,0)) to lo.bound
  start_area = tryCatch(
    pft(upper = max(c(lo.bound - 10*sigma,0)), df = df, mu = abs(mu), 
        sigma = sigma, lower = 0) +
      pft(upper = lo.bound, df = df, mu = abs(mu), 
          sigma = sigma, lower = max(c(lo.bound - 10*sigma,0))),
    error = function(c) NaN)
  
  # Then we integrate from xbar_dm and above to for the root, added the start_area
  # to the integration
  xroot <-  tryCatch(
    uniroot( function(z)
      pft(upper = z, df = df, mu = abs(mu), sigma = sigma, lower = lo.bound) + 
        start_area - p,
      lower = lo.bound, upper = hi.bound, extendInt = "no")$root,
    error = function(c) NaN)

  # browser();
  
  return(xroot)
}



mdm_tdist <- function(x, y = NULL, conf.level = 0.95) {
  #' @description Calculate most difference in means statistic from integrating 
  #' over a folded t-distribution to an extent dicated by the conf.level
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @return Returns most difference in means (single value)
  #' @usage mdm_tdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=3,mean=0,sd=1); y <- rnorm(n=3,mean=0,sd=1);
  #' mdm_tdist(x,y, conf.level = 0.95)

  # save(list = ls(all.names = TRUE), file = "temp/mdm_tdist.RData",envir = environment())
  # load(file = "temp/mdm_tdist.RData")

  # Calculate basic stats x and y
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)

  # Calculate difference in means stats
  if (is.null(y)) {
    # 1-sample case
    # Degrees of freedom and upper search bound for quantile function calculated 
    # from built-in welch's t-test
    wtt <- t.test(x=x,y=NULL, var.equal = FALSE, conf.level = 1-(1-conf.level)/4)
    df_d <- wtt$parameter
    
    xbar_dm <- mean(x)
    sd_d <- sd_x
    sd_dm = sd_x / sqrt(n_x)

  } else {
    # 2-sample case
    # Degrees of freedom and upper search bound for quantiel function calculated 
    # from built-in welch's t-test
    wtt <- t.test(x=y,y=x, var.equal = FALSE, conf.level = 1-(1-conf.level)/4)
    df_d <- wtt$parameter

    xbar_dm <- mean(y) - mean(x)
    sd_d <- sqrt( (( n_x - 1) * sd_x^2  +  (n_y - 1) * sd_y^2 ) / df_d)
    sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
  }

  # Calculate search bounds for mdm, used for uniroot call
  lo.bound= abs(xbar_dm)
  hi.bound = max(abs(wtt$conf.int))
  if (lo.bound>hi.bound) {browser();}

  # Quantile with prob set to alpha and sampple estimates as parameters
  mdm <-  tryCatch(qft(p = conf.level, df = df_d, mu = xbar_dm, sigma = sd_dm,
                       lo.bound, hi.bound), error = function(c) NaN)
  
  # If integration fails pause execution
  if (is.nan(mdm) || is.null(mdm)) {browser();}
  
  return(mdm)
}


# 
# df_WelchSatter<- function(x,y) {
#   s1 = sd(x)
#   s2 = sd(y)
#   n1 = length(x)
#   n2 = length(y)
# 
#   v1 = n1 - 1; v2 = n2 - 1
#   v = (s1^2/n1 + s2^2/n2)^2/ ( s1^4/(n1^2/v1)  +  s2^4/(n2^2/v2)  )
#   return(v)
# }




macb_tdist_2sample <- function (x, y, conf.level = 0.95) {
  #' Most absolute one-tailed confidence bounds of t-distribution
  #' 
  lo_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "greater")$conf.int[1]
  up_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "less")$conf.int[2]
  
  macb <- max(abs(c(lo_b, up_b)))
  
  return(macb)
}


lacb_tdist_2sample <- function (x, y, conf.level = 0.95) {
  #' Least absolute one-tailed confidence bounds of t-distribution
  #' 
  lo_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "greater")$conf.int[1]
  up_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "less")$conf.int[2]

  lacb <- min(abs(c(lo_b, up_b)))
  if (sign(ci[1])!=signci[2]) {lacb <- 0}
  
  return(lacb)
}



  
rmdm_tdist <- function(x_ctrl, y_exp, conf.level = 0.95, 
                              verbose = FALSE,  var.equal = FALSE, method = "fieller")  {
  #' @description Calculates the relative most difference in means assuming with
  #' rmdm = mdm/X, X being the control group and Y the experimental
  #' 
  #' @param x_ctrl vector of measurements in control group
  #' @param y_exp vector of measurements in experimental group
  #' @param conf.level significance level for calculating upper mdm
  #' 
  #' @return relative most difference in means
  
  # Equation for pooled variance taken from:
  # https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals5.html

  # Calculate basic sample statistics
  mean_ctrl = mean(x_ctrl)
  sd_ctrl = sd(x_ctrl)
  n_ctrl = length(x_ctrl)
  se_ctrl = sd_ctrl/sqrt(n_ctrl)
  
  mean_exp = mean(y_exp)
  sd_exp = sd(y_exp)
  n_exp = length(y_exp)
  se_exp = sd_exp/sqrt(n_exp)
  
  mean_dm = mean_exp - mean_ctrl
  
  # sd_dm is equivalent to se_d (standard error of difference)
  if (var.equal) {
    sd_dm = sqrt(  ( (n_ctrl - 1)*sd_ctrl^2 + (n_exp - 1)*sd_exp^2 ) / (n_ctrl + n_exp - 2)  ) * sqrt(1/n_exp + 1/n_exp)
  } else {
    sd_dm = sqrt(sd_ctrl^2/n_ctrl + sd_exp^2/n_exp)
  }

  # Calculate RMDM based on specified method
  if (method =="qnormrat") {
    # # Each one tailed confidence bound specified by conf.level
    # qn_lo <- qnormrat(1-conf.level, mean_dm, sd_dm, abs(mean_ctrl), se_ctrl)
    # qn_hi <- qnormrat(conf.level, mean_dm, sd_dm, abs(mean_ctrl), se_ctrl)
    # rmdm = max(abs(c(qn_lo, qn_hi)))
    rmdm <- qnormrat(conf.level, abs(mean_dm), sd_dm, abs(mean_ctrl), se_ctrl)
    
    
  } else if (method =="fieller") {

    # browser()
    f <- n_ctrl + n_exp - 2
    sSquared <- (sum((x_ctrl - mean_ctrl)^2) + sum((y_exp - mean_exp)^2))/f
    
    fieller_int <- tryCatch(get_FiellerInterval(mean_ctrl, mean_exp, sSquared, 
                                                1/n_ctrl, 1/n_exp, f, v12 = 0, alpha=1-conf.level), 
                            error = function(c) data.frame(upper = NaN, lower = NaN))
    
    rmdm <- max(abs(c(fieller_int$lower,fieller_int$upper)))

  } else {stop('rmdm_tdist(): unsupported method')}
  
    
  # browser()
  
  return(rmdm)
}


