
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


mdm_normal <- function(x, y = NULL, paired = FALSE, var.equal = FALSE, conf.level = 0.95, 
           verbose = FALSE, distribution = NULL) {
  #' Calculate most difference in means using a normal distribution (z or t)
  #'
  #' @description  most difference in means is the largest 
  #' difference than could exist between group(s) in the data.
  #' 
  #' Calculates most difference in means assuming a normal distribution
  #' (z or t) by integrating a central normal pdf shifted by sample mean and
  #' Using a root finding function. most difference in means is the largest 
  #' difference than could exist between group(s) in the data.
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param paired boolean that specifies if 2-sample data is paired (default: FALSE)
  #' @param var.equal boolean if variance assumed equal (default: FALSE)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @param verbose function prints out extra input during execution (default: FALSE)
  #' @param distribution string specifying assumed distribution ('t_dist', 
  #' 'z_dist') (Default: NULL). If not specified, t_dist chosen for cases under 
  #' 30 samples
  #' @return Returns most difference in means (single value)
  #' @usage mdm_normal_tdist(x)
  #' @usage mdm_normal_tdist(x, y)
  #' @usage mdm_normal_tdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=10,mean=0,sd=1); 
  #' y <- rnorm(n=10,mean=0,sd=1); 
  #' mdm_normal(x,y)
  # Quick error check 
  if (is.null(x)) errorCondition('input argument x must be specified.')
  
  # Convert 2-sample paired to 1-sample
  if (paired & !is.null(y)) {x = x-y; y=NULL}
  
  # Calculate degrees of freedom and 
  if (is.null(y)) df <- length(x)-1 else df <- length(x) + length(y) - 2
  # Determine Distribution if not specified, n=30 threshold between z and t dist.
  if (is.null(distribution) & df<=30) {
    distribution <- 't_dist' } else {distribution<-'z_dist'}
  if (verbose) print(sprintf('Distribution: %s', distribution))
  
  # Calculate MDM with specified distribution
  if (distribution=='t_dist'){
    mdm <- mdm_normal_tdist(x, y, conf.level = conf.level, verbose = verbose,
                                     search_pad_percent = 0.1,var.equal = var.equal)
  } else if (distribution =='z_dist') {
    mdm = mdm_normal_zdist(x,y, conf.level = conf.level, verbose = verbose,
                                search_pad_percent = 0.1,var.equal = var.equal)
  } else {
    mdm=NULL
    warning("Input argument 'distribution' not supported")
  }
  
  return(mdm)
}







dft <- function(x, df, mu, sigma) {
  #' @description: probability density function

  # save(list = ls(all.names = TRUE), file = "temp/dft.RData",envir = environment())
  # load(file = "temp/plot_population_params.RData")
  
  v = df
  t1 <- (1 +  (1/v) * ((x - mu)^2 / sigma^2 )) ^ ( -(v+1) / 2)
  t2 <- (1 +  (1/v) * ((x + mu)^2 / sigma^2 )) ^ ( -(v+1) / 2)
  
  d <- gamma( (v+1)/2 ) / (gamma(v/2) * sqrt(v*pi*sigma^2) ) * (t1 + t2)
  # d is zero for all x<0 by definition
  d[x<0]<- 0

  return(d)
}

pft <- function(upper, df, mu, sigma, lower = -Inf) {
  #' @description: cumulative distribution function
  # save(list = ls(all.names = TRUE), file = "temp/pft.RData",envir = environment())
  # load(file = "temp/plot_population_params.RData")
  
  # browser();
  
  # start_time <- Sys.time()
  p <- integrate(f = function(y) dft(y, df, abs(mu),sigma),
                 lower = lower, upper = upper,
                  rel.tol = 1e-8, abs.tol = 1e-9)$value
  # end_time <- Sys.time()
  # end_time - start_time
  # 
  # start_time <- Sys.time()
  # p <- cubintegrate(f = function(y) dft(y, df, abs(mu),sigma),
  #                lower = lower, upper = upper,
  #                relTol = 1e-6, absTol = 1e-12)$integral
  # end_time <- Sys.time()
  # end_time - start_time
  
  
  return(p)
}

qft <- function(p, df, mu, sigma) {
  #' @description Quantile function of ratio normal distribution (X/Y)
  #' For given probability, return position of that percentile from cdf
  #' 
  #' @param p probability
  # Estimate bounds to search for root
  
  # save(list = ls(all.names = TRUE), file = "temp/qft.RData",envir = environment())
  # load(file = "temp/plot_population_params.RData")
  
  lo.bound = abs(mu)
  hi.bound = abs(mu) + qt(1-(1-p)/4, df) * sigma 
  
  # x = seq(lo.bound, hi.bound, (hi.bound - lo.bound)/1e4)
  # gg <- ggplot(data = data.frame(x=x,y=sapply(x,function(x)
  #   pft(x, df = 100, mu = mu, sigma = sigma))),
  #              aes(x=x, y=y)) + geom_line()
  # print(gg)
  
  start_area = pft(x=lo.bound, df=df, mu=abs(mu), sigma=sigma, start = 0)
  
  xroot <-  uniroot(function(z) 
    pft(x=z, df=df, mu=abs(mu), sigma=sigma, start = lo.bound) + start_area -p, 
    lower = lo.bound, upper = hi.bound, extendInt = "no", tol = .Machine$double.eps^.5)
  
  return(xroot$root)
}








mdm_normal_zdist <- function(x, y = NULL, conf.level = 0.95, verbose = FALSE, 
                             var.equal = FALSE, search_pad_percent = 0.01, 
                             method="nonstandard") {
  #' Calculate most difference in means using z distribution
  #'
  #' @description Calculate most difference in means statistic from integrating a 
  #' central normal pdf shifted by -x_bar Using root finding function to 
  #' integrate over a normal CDF with area under the curve equal to (1-a). 
  #' Calculated from the difference in means distribution for two samples, 
  #' or keeps same distribution for one sample.
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @param verbose function prints out extra input during execution
  #' @param var.equal boolean assume variance qual between groups (TODO: not 
  #' implemented for TRUE)
  #' @param search_pad_percent for calculating statistic, specifies how outside 
  #' the theoretical search window the root finding can look.
  #' @return Returns most difference in means (single value)
  #' @usage mdm_normal_zdist(x)
  #' @usage mdm_normal_zdist(x, y)
  #' @usage mdm_normal_zdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=50,mean=0,sd=1); y <- rnorm(n=50,mean=0,sd=1); 
  #' mdm_normal_zdist(x,y)
  
  # save(list = ls(all.names = TRUE), file = "temp/mdm_normal_zdist.RData",envir = environment())
  # load(file = "temp/mdm_normal_zdist.RData")
  
  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distirbution of the sampel for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample stats
    wtt <- t.test(x=x,y=NULL, var.equal = FALSE, conf.level = 1-(1-conf.level)/4)
    df_d <- wtt$parameter
    xbar_dm <- mean(x)
    sd_d <- sd_x
    sd_dm = sd_x / sqrt(n_x)
    
  } else { 
    # 2-sample stats
    # df_d <- df_WelchSattern(x,y)
    wtt <- t.test(x=y,y=x, var.equal = FALSE, conf.level = 1-(1-conf.level)/4)
    df_d <- wtt$parameter
    
    xbar_dm <- mean(y) - mean(x)
    sd_d <- sqrt( (( n_x - 1) * sd_x^2  +  (n_y - 1) * sd_y^2 ) / df_d)
    sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
    # Assume unequal variance
  }
  

  # browser();
  
  lo.bound= abs(xbar_dm)
  hi.bound = max(abs(wtt$conf.int))
  if (lo.bound>hi.bound) {browser();}
  
  
  
  # The MDM is always greater than the sample mean, so integrate to that point first
  # Integration is broken up bceause in cases where xbar_dm >> sd_dm, the pdf is 
  # so sparse that integration algorithm can miss the peak
  start_area = tryCatch(
    pft(upper = max(c(lo.bound - 10*sd_dm,0)), df = df_d, mu = abs(xbar_dm), sigma = sd_dm, lower = 0) +
    pft(upper = lo.bound, df = df_d, mu = abs(xbar_dm), sigma = sd_dm, max(c(lo.bound - 10*sd_dm,0))), 
    error = function(c) NaN)
  
  mdm <-  tryCatch(
    uniroot( function(z) 
    pft(upper=z, df=df_d, mu=abs(xbar_dm), sigma=sd_dm, lower = lo.bound) + start_area - conf.level, 
    lower = lo.bound, upper = hi.bound, extendInt = "no")$root, 
                            error = function(c) NaN)
  
  #todo: add mu.sigma cuttoff where abs value transform ignored (will get greater accuracy)


  # x = seq(lo.bound, hi.bound, (hi.bound - lo.bound)/1e4)
  # gg <- ggplot(data = data.frame(x=x,y=sapply(x,function(x)
  #   pft(x, df = 100, mu = abs(xbar_dm), sigma = sd_dm))),
  #              aes(x=x, y=y)) + geom_line()
  # print(gg)
  # 
  # 
  # x = seq(lo.bound-0.5, hi.bound+2, (hi.bound - lo.bound)/1e4)
  # gg <- ggplot(data = data.frame(x=x,y=sapply(x,function(x)
  #   dft(x, df = 100, mu = abs(xbar_dm), sigma = sd_dm))),
  #   aes(x=x, y=y)) + geom_line()
  # print(gg)


  # mdm <- tryCatch(qft(p = conf.level, df = df_d, mu = xbar_dm, sigma = sd_dm), 
  #                         error = function(c) NaN)
  if (is.nan(mdm) || is.nan(start_area)) {browser();}
  
  # mdm = qft(conf.level, df_d, abs(xbar_dm), sd_dm)
  
  return(mdm)  
}


df_WelchSatter<- function(x,y) {
  s1 = sd(x)
  s2 = sd(y)
  n1 = length(x)
  n2 = length(y)
  
  v1 = n1 - 1; v2 = n2 - 1
  v = (s1^2/n1 + s2^2/n2)^2/ ( s1^4/(n1^2/v1)  +  s2^4/(n2^2/v2)  )
  return(v)
}

check_bounds <- function(x, search_bounds, verbose = FALSE, range_tol=1000) {
  
  is_warn=FALSE;
  if ( x > search_bounds[2] + abs(diff(search_bounds)/range_tol) ) {
    warning("mdm: equation root equal to upper bounds of search space: 
            results unreliable."); is_warn = TRUE;
  }
  
  if (x < search_bounds[1] - abs(diff(search_bounds)/range_tol)   ) {
    warning("mdm: equation root equal to lower bounds of search space: 
            results unreliable."); is_warn = TRUE;
  }
  
  if (is_warn) {browser();}
  return(is_warn)
}


ucl_tdist_mean <- function(x_bar, sx, n_x, dfx, semx = NULL, alpha = 0.05) {
  #  Calculate t statistic from integrating a central normal pdf shifted by -x_bar
  # Using minimization function to integrate using a normal CDF setting the area 
  # under the curve to (1-a)
  # Inputs
  #   x_bar: sample mean
  #   sx: Sample standard deviation
  #   n_x: sample number
  #   dfx: sample degrees of freedom
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Estimate multiplier for bounds of search interval
  sd_mult <- qt(1 - alpha, dfx)
  
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  t_star_function <- function(x) { abs( (pt(x,df = dfx) - pt(-x, df = dfx))  - (1 - alpha))}
  
  # Calculate roughly the extext of search space
  search_bounds = c(0, x_bar + 4*sd_mult * sx)
  
  # Minimization
  t_star = optimize(t_star_function, search_bounds, maximum = FALSE)
  
  if (t_star$minimum == search_bounds[2]) {
    warning("mdm: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = s_x / sqrt(n_x)}
  
  # CI_mean x_bar +- t_star * s/sqrt(n)
  ucl = x_bar + t_star$minimum * sem_x
  
  print(t_star$minimum)
  
  return(ucl)
}



max_abs_cl_mean_z <-  function(x_bar, sem_x, alpha, x=NULL) {
  
  if (!is.null(x)) {
    x_bar <- mean(x); sem_x <- sd(x)/sqrt(length(x))
  }
  
  conf.lims <- qnorm( c(alpha/2, 1 - alpha/2), mean = x_bar, sd = sem_x)
  max_abs_cl <- max(abs(conf.lims))
  return(max_abs_cl)
}


max_abs_cl_mean_z_standard <-  function(x_bar, sem_x, alpha) {
  cl_mean <- qnorm( c(alpha/2, 1 - alpha/2), sd = sem_x)
  max_abs_cl <- max(abs(cl_mean))
  return(max_abs_cl)
}






# Need better method to calculate ratio 95% ci
# https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/ratimean.htm
# https://www.graphpad.com/quickcalcs/errorProp2/
# https://i.stack.imgur.com/vO8Ip.png


# Discussion
# https://stats.stackexchange.com/questions/16349/how-to-compute-the-confidence-interval-of-the-ratio-of-two-normal-means
# https://www.rdocumentation.org/packages/mratios/versions/1.3.17/topics/t.test.ratio

# Relative form of mdm
  
rmdm_normal_zdist <- function(x_ctrl, y_exp, conf.level = 0.95, 
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

  } else {stop('rmdm_normal_zdist(): unsupported method')}
  
    
  # browser()
  
  return(rmdm)
}



# rmdm_normal_montecarlo <- 
#   function(m_x, se_x, m_y, se_y, means_x = NULL, means_y = NULL, conf.level = 0.95,
#            n_trials=1e6)  {
# 
#   
#   if (m_y < 0) {stop('Mean of y must be > 0')}
#   
#   # Each trials just draws means instead of generating data points
#   if (is.null(means_x)){ means_x = rnorm(n_trials, mean = m_x, sd = se_x)}
#   if (is.null(means_y)){ means_y = rnorm(n_trials, mean = m_y, sd = se_y)}
#   
#   
#   stat = abs(means_y-means_x)/means_y
#   
#   rmdm = quantile(stat, conf.level)
#   
# }

