
# MDM Package: most difference in means
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.


# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(docstring)
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


mdm_normal_tdist <- function(x,y = NULL, conf.level = 0.95, verbose = FALSE, 
                             var.equal = FALSE, search_pad_percent = 0.01) {
  #' Calculate most difference in means using t distribution
  #'
  #' @description  Calculate most difference in means statistic from integrating a
  #' central normal t-distribution pdf shifted by -x_bar Using a root finding 
  #' function to integrate over a normal CDF, with area under the curve equal 
  #' to (1-a). Calculate stats of the difference in means distribution for two 
  #' samples, or keeps same distribution for one sample case.
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @param verbose function prints out extra input during execution
  #' @param search_pad_percent for calculating statistic, specifies how far 
  #' outside the theoretical search window the root finding can look.
  #' @return Returns most difference in means (single value)
  #' @usage mdm_normal_tdist(x)
  #' @usage mdm_normal_tdist(x, y)
  #' @usage mdm_normal_tdist(x, y, conf.level = 0.95)
  #' @examples 
  #' x <- rnorm(n=8,mean=0,sd=1); 
  #' y <- rnorm(n=8,mean=0,sd=1); 
  #' mdm_normal_tdist(x,y)

  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distirbution of the sampel for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample stats
    xbar_dm <- mean(x)
    df_d <- n_x - 1
    sd_d <- sd_x
    sd_dm = sqrt( sd_x^2 / n_x)
    
  } else { 
    # 2-sample stats
    xbar_dm <- mean(y) - mean(x)
    df_d <- n_x + n_y - 2
    sd_d <- sqrt( (( n_x-1)*sd_x^2 + (n_y - 1) * sd_y^2 ) / df_d)
    sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
    if (verbose) print(sprintf("xbar_dm: %.3f", xbar_dm))  
  }
  
  # Calculate search bounds defined by tails of alpha and 2*alpha CI of mean 
  #alpha = (1 - conf.level)
  #ci_mean_alpha  <- qt( c(  alpha/2, 1 -   alpha/2), ncp = xbar_dm, sd = sd_d)
  #ci_mean_2alpha <- qt( c(2*alpha/2, 1 - 2*alpha/2), ncp = xbar_dm, sd = sd_d)
  
  # Calculate search bounds defined by tails of alpha and 2*alpha CI of mean 
  ci_mean_alpha  <- t.test(x,y, conf.level = conf.level, paired = FALSE, 
                       var.equal = FALSE, alternative = "two.sided")$conf.int
  ci_mean_2alpha <- t.test(x,y, conf.level = 1-2*(1-conf.level), paired = FALSE, 
                       var.equal = FALSE, alternative = "two.sided")$conf.int
  lower_bounds =  max(abs(ci_mean_2alpha))
  upper_bounds = max(abs(ci_mean_alpha))
  # Add extra padding around search bounds so root finding not done on boundary
  bounds_range = upper_bounds - lower_bounds
  search_bounds = c(lower_bounds - search_pad_percent * bounds_range,
                    upper_bounds + search_pad_percent * bounds_range)
  if (verbose) print(sprintf('Bounds:[ %.3f  %.3f]', search_bounds[1], search_bounds[2]))
  
  # Calculate MDM with root finding optmization
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  t_star_standard <- function(x) {pt(q = (-xbar_dm + x) / sd_dm, df = df_d) - 
                                  pt(q = (-xbar_dm - x) / sd_dm, df = df_d) - conf.level}
  
  # Debugging search bounds are correct
  if (verbose) {
    t = seq(from = search_bounds[1], to = search_bounds[2], by = diff(search_bounds)/100)
    f_t = sapply(t, t_star_standard)
    plot(t,f_t)
  }
  
  
  # Solve for t star with root finding where t_star_standard equals (1 - alpha)
  t_star = uniroot(t_star_standard, search_bounds, check.conv = TRUE,
                         tol = .Machine$double.eps, maxiter = 1000, trace = 0)
  # The optimized root should fall entirely within the earch bounds 
  check_bounds(t_star, search_bounds, verbose = verbose, range_tol = 1000)
    
  # MDM is root location added to difference of means
  if (verbose) print(sprintf("t_star: %.3f", t_star$root))
  mdm = t_star$root 
  
  return(mdm)
}



mdm_normal_zdist_approx <- function(x, y = NULL, conf.level = 0.95) {
  #' Calculate most difference in means using folded z distribution
  #'
  #' @description Calculate most difference in means statistic from a folded 
  #' normal z distribution. Uses R built-in qfoldnorm to calculate, but this 
  #' function is an approximation that has noticeably worse coverage error 
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
  #' @usage mdm_normal_zdist_approx(x)
  #' @usage mdm_normal_zdist_approx(x, y)
  #' @usage mdm_normal_zdist_approx(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=50,mean=0,sd=1); y <- rnorm(n=50,mean=0,sd=1); 
  #' mdm_normal_zdist_approx(x,y)
  
  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distirbution of the sampel for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample stats
    df_d <- n_x - 1
    xbar_dm <- mean(x)
    sd_d <- sd_x
    sd_dm = sd_x / sqrt(n_x)
    
  } else { 
    # 2-sample stats
    df_d <- n_x + n_y - 2
    xbar_dm <- mean(y) - mean(x)
    sd_d <- sqrt( (( n_x - 1) * sd_x^2  +  (n_y - 1) * sd_y^2 ) / df_d)
    sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
    # if (verbose) print(sprintf("xbar_dm: %.3f", xbar_dm))  
  }

  # Calculate the quantile function for folded normal with built-in qfoldnorm
  mdm <- qfoldnorm(conf.level, mean = xbar_dm, sd = sd_dm, a1 = 1, a2 = 1,
            lower.tail = TRUE, log.p = FALSE)

  return(mdm)  
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
  
  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distirbution of the sampel for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample stats
    df_d <- n_x - 1
    xbar_dm <- mean(x)
    sd_d <- sd_x
    sd_dm = sd_x / sqrt(n_x)
    
  } else { 
    # 2-sample stats
    df_d <- n_x + n_y - 2
    xbar_dm <- mean(y) - mean(x)
    sd_d <- sqrt( (( n_x - 1) * sd_x^2  +  (n_y - 1) * sd_y^2 ) / df_d)
    sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
    # if (verbose) print(sprintf("xbar_dm: %.3f", xbar_dm))  
  }
  
  # Calculate search bounds defined by tails of alpha and 2*alpha CI of mean 
  alpha = (1 - conf.level)
  if (method =="standard") {
    lower_bounds = max_abs_cl_mean_z_standard(xbar_dm, sd_dm, 2*alpha)
    upper_bounds= max_abs_cl_mean_z_standard(xbar_dm, sd_dm, alpha)
    
    # Note: xbar_dm and x do not need to be in z_score units because they are z 
    # normalized after insertion into function
    z_star_fcn <- function(x) { pnorm( (-xbar_dm + x)/sd_dm, mean = 0, sd = 1) - 
        pnorm( (-xbar_dm - x)/sd_dm, mean = 0, sd = 1) - conf.level - xbar_dm}
    
  } else if (method =="nonstandard") {
    lower_bounds = max_abs_cl_mean_z(xbar_dm, sd_dm, 2*alpha)
    upper_bounds = max_abs_cl_mean_z(xbar_dm, sd_dm, alpha)
    
    # Integration of folded z-distribution from standard central z-distribution
    z_star_fcn <- function(x) {pnorm( +x, mean = xbar_dm, sd = sd_dm) - 
        pnorm( -x, mean = xbar_dm, sd = sd_dm) - conf.level}
    
  } 

  # Add extra padding around search bounds for root finding at boundary
  bounds_range = upper_bounds - lower_bounds
  
  
  search_bounds = c(lower_bounds - search_pad_percent * bounds_range,
                    upper_bounds + search_pad_percent * bounds_range)
  # if (verbose) print(sprintf('Bounds:[ %.3f  %.3f]', search_bounds[1], search_bounds[2]))

  # Solve for t star with root finding
  z_star = uniroot(z_star_fcn, search_bounds, check.conv = TRUE,
                   tol = .Machine$double.eps, maxiter = 1000, trace = 0, 
                   extendInt="upX")
  
  if (verbose) {
    z = seq(from = search_bounds[1], to = search_bounds[2], by = diff(search_bounds)/100)
    f_z = lapply(z, z_star_fcn)
    plot(z,f_z,type="l"); abline(v=z_star$root,col="blue")
    abline(v = lower_bounds, col="red", lwd=3, lty=2)
    abline(v = upper_bounds, col="red", lwd=3, lty=2)
    z_star$root-lower_bounds
  }
  
  
  # The optimized root should fall entirely within the earch bounds 
  is_warn = check_bounds(z_star$root, c(lower_bounds, upper_bounds), verbose = FALSE, range_tol= 1000)
  if (is_warn) {browser()}
  
  # MDM is root of integration
  if (verbose) print(sprintf("z_star: %.4e, lower:%.4e, upper:%.4e", z_star$root,lower_bounds,upper_bounds))
  mdm = z_star$root
  

  return(mdm)  
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
    
    fieller_int <- get_FiellerInterval(mean_ctrl, mean_exp, sSquared, 1/n_ctrl, 1/n_exp, f, v12 = 0, alpha=1-conf.level)
    rmdm <- max(abs(c(fieller_int$lower,fieller_int$upper)))

  } else {stop('rmdm_normal_zdist(): unsupported method')}
  
    
  # browser()
  
  return(rmdm)
}



rmdm_normal_montecarlo <- 
  function(m_x, se_x, m_y, se_y, means_x = NULL, means_y = NULL, conf.level = 0.95,
           n_trials=1e6)  {

  
  if (m_y < 0) {stop('Mean of y must be > 0')}
  
  # Each trials just draws means instead of generating data points
  if (is.null(means_x)){ means_x = rnorm(n_trials, mean = m_x, sd = se_x)}
  if (is.null(means_y)){ means_y = rnorm(n_trials, mean = m_y, sd = se_y)}
  
  
  stat = abs(means_y-means_x)/means_y
  
  rmdm = quantile(stat, conf.level)
  
}

