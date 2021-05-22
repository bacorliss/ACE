
# MDM Package: most difference in means
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.


# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(docstring)




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

# max_abs_t_cl_mean <-  function(x, df, alpha) {
#   print("max_abs_t_cl_mean: Warning: Untested")
#   cl_mean <- mean(x) + c(
#     qnorm(1-(alpha/2), df = length(x)-1) * sd(x)/sqrt(length(x)),
#     qnorm(   alpha/2,  df = length(x)-1) * sd(x)/sqrt(length(x)))
#   max_abs_cl = max(abs(cl_mean))
# }



# Utility functions for ratio normal distribution
# dnormrat: density function
# qnormrat: cumul density function
# qnormrat: quantile function

# Error function
erf <- Vectorize(function(x) 2*pnorm(sqrt(2)*x) - 1)

dnormrat <- function(z, mu_x, sigma_x,mu_y, sigma_y) {
  #' Probability density function of ratio of two normal distributions
  # Based on code from: 
  # https://rstudio-pubs-static.s3.amazonaws.com/287838_7c982110ffe44d1eb5184739c5724926.html
  
  
  erf <- Vectorize(function(x) 2*pnorm(sqrt(2)*x) - 1)
  
  beta = mu_x/mu_y
  rho = sigma_y/sigma_x
  deltay =  sigma_y/mu_y
  
  q <- (1+beta*rho^2*z)/(deltay*sqrt(1+rho^2*z^2))
  fz <- rho/(pi*(1+rho^2*z^2))*( exp(-(rho^2*beta^2+1)/(2*deltay^2))  + 
                                    sqrt(pi/2)*q*erf(q/sqrt(2))*exp(-0.5*(rho^2*(z-beta)^2)/(deltay^2*(1+rho^2*z^2))) )
   return(fz)  
}

pnormrat <- function(z, mu_x, sigma_x,mu_y, sigma_y, start = -Inf) {
  #' Cumulative probability density function
  #' 
  Fz <- integrate(function(x) dnormrat(x, mu_x, sigma_x,mu_y, sigma_y), start, z, 
                  rel.tol = .Machine$double.eps^.5, abs.tol=.Machine$double.eps^.5)
  
  return(Fz$value)
}

qnormrat <- function(p, mu_x, sigma_x, mu_y, sigma_y, VERBOSE=FALSE) {
  # For given probability, find position of that percentile from cdf

 # Calculate conservative search bounds for finding the 95 percentile area
 lo.bound = (mu_x - qnorm(p)*sigma_x) / (mu_y + qnorm(p)*sigma_y)
 hi.bound = (mu_x + qnorm(p)*sigma_x) / (mu_y - qnorm(p)*sigma_y)
 
 # To speed up uniroot compuation, we are not going to integrate from -Inf with 
 # every call to pnormrat, instead, we will calculate lower bound area and restrict 
 # our uniroot integration from the area of the lower search bound. So pnormrat 
 # does not have to integrate from -Inf with every function call.
 
 # Calculate area to lo.bound
 lo.bound.p = pnormrat(lo.bound, mu_x, sigma_x, mu_y, sigma_y)
 
  qroot <- uniroot(function(z) 
    pnormrat(z, mu_x, sigma_x, mu_y, sigma_y, start = lo.bound) + lo.bound.p - p, 
          lower = abs(lo.bound), upper = abs(hi.bound), extendInt = "no")
  if (VERBOSE) {
    sprintf("root:%f [%f %f]",qroot$root, lo.bound, hi.bound)
  }
  # 
  
  # if (qroot$root<lo.bound || qroot > hi.bound) {print("root outside of search zone")
  #   browser();}
  return(qroot$root)
}

# qnormrat2 <- function(p, mu_x, sigma_x, mu_y, sigma_y, VERBOSE=FALSE) {
#   # For given probability, find position of that percentile from cdf
#   
#   lo.bound = (mu_x - qnorm(p)*sigma_x) / (mu_y + qnorm(p)*sigma_y)
#   hi.bound = (mu_x + qnorm(p)*sigma_x) / (mu_y - qnorm(p)*sigma_y)
#   
#   # Calculate area to lo.bound
#   # lo.bound.p = pnormrat(z, mu_x, sigma_x, mu_y, sigma_y)$value
#   
#   # Calculate area 
#   
#   # Estimate bounds
#   qroot <- uniroot(function(z) pnormrat(z, mu_x, sigma_x, mu_y, sigma_y) - p, 
#                    lower = abs(lo.bound), upper = abs(hi.bound), extendInt = "no")
#   if (VERBOSE) {
#     sprintf("root:%f [%f %f]",qroot$root, lo.bound, hi.bound)
#   }
#   # 
#   
#   # if (qroot$root<lo.bound || qroot > hi.bound) {print("root outside of search zone")
#   #   browser();}
#   return(qroot$root)
# }



# Need better method to calculate ratio 95% ci
# https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/ratimean.htm
# https://www.graphpad.com/quickcalcs/errorProp2/
# https://i.stack.imgur.com/vO8Ip.png


# Discussion
# https://stats.stackexchange.com/questions/16349/how-to-compute-the-confidence-interval-of-the-ratio-of-two-normal-means
# https://www.rdocumentation.org/packages/mratios/versions/1.3.17/topics/t.test.ratio

# Relative form of mdm
  
rmdm_normal_zdist <- function(x, y, mdm = NULL, conf.level.mdm = 0.95, conf.level.rmdm = 0.95, 
                              verbose = FALSE,  var.equal = FALSE)  {
  #' @description Calculates the relative most difference in means assuming with
  #' rmdm = mdm/X, X being the control group and Y the experimental
  #' 
  #' @param x vector of measurements in control group
  #' @param y vector of measurements in experimental group
  #' @param conf.level.mdm significance level for calculating upper mdm
  #' @param conf.level.rmdm significance level for calculating upper rmdm
  #' 
  #' @return relative most difference in means
  

  # Calculate mean and std of difference in means distribution
  xbar_dm = mean(y) - mean(x) 
  s_dm <- sqrt(var(x)/length(x) + var(y)/length(y)) 
  
  # Calculate new mean and std after absolute value transform
  df <- moments_foldnorm(xbar_dm, s_dm)
  
  # Calculate upper percentile specified by 
  # rmdm <- qnormrat(conf.level.rmdm, df$mu_f, df$sigma_f, mean(x), sd(x)/sqrt(length(x)))
  
  
  
  
  # Calculate 95% CI of ratio of means, must input std into equations
  # rmdm <- ttest_diffratio(df$mu_f, sqrt(sd(y)^2 + sd(x)^2), length(x)+length(y)-1, 
  #                         mean(x), sd(x), length(y), alternative = "less", rho = 1, 
  #                              var.equal = FALSE, conf.level = 0.95)$conf.int[2]
    
    
  # browser()
  
  return(rmdm)
}


moments_foldnorm <- function(mu, sigma) {
  #' @description given the mean and standard deviation of a normal distribution 
  #' before transform, calculates the mean, std, and variance after absolute 
  #' value transform.
  #' 
  #'  Equations taken from
  #'    https://www.randomservices.org/random/special/FoldedNormal.html
  #' Tested against Special Distribution Simulator
  #'    https://www.randomservices.org/random/apps/SpecialSimulation.html
  
  # Equations tested against
  mu_f = mu * (1 - 2*pnorm(-mu/sigma)) + sigma*sqrt(2/pi)* exp(-mu^2/(2*sigma^2))
  # Ex2 = mu^2 + sigma^2
  var_f = mu^2 + sigma^2  - mu_f^2
  
  
  return(data.frame(mu_f=mu_f,sigma_f = sqrt(var_f), var_f = var_f))
}





# 
# 
# ttest_diffratio <- function (mean_x, sd_x, nx, mean_y, sd_y, ny, alternative = "two.sided", rho = 1, 
#                              var.equal = FALSE, conf.level = 0.95, iterativeCI=FALSE,
#                              ul=1e10, ll=-1e10, ...) 
# { # code modified from:
#   #https://cran.r-project.org/web/packages/mratios/index.html
#   
#   addargs <- list(...)
#   alternative <- match.arg(alternative, choices = c("two.sided", 
#                                                     "less", "greater"))
#   if (!is.numeric(rho) | length(rho) != 1) {
#     stop("Argument 'rho' must be a single numeric value")
#   }
#   if (!is.logical(var.equal) | length(var.equal) != 1) {
#     stop("Argument'var.equal' must be either TRUE or FALSE")
#   }
#   if (!is.numeric(conf.level) | length(conf.level) != 1 | conf.level <= 
#       0.5 | conf.level >= 1) {
#     stop("Argument 'conf.level' must be a single numeric value between 0.5 and 1")
#   }
# 
#   mx <- mean_x
#   my <- mean_y
#   vx <- sd_x^2
#   vy <- sd_y^2
#   est <- mx/my
#   
#   if (sqrt(vx) < 10 * .Machine$double.eps * abs(mx)) {
#     stop("data in x are essentially constant")
#   }
#   if (sqrt(vy) < 10 * .Machine$double.eps * abs(my)) {
#     stop("data in y are essentially constant")
#   }
#   if (is.null(addargs$namex) || is.null(addargs$namey)) {
#     namex = "x"
#     namey = "y"
#   }
#   else {
#     namex = addargs$namex
#     namey = addargs$namey
#   }
#   if(any(c(my,mx)<0)){warning("Sample means are smaller than 0! References for this test do not consider this case explicitly!")}
#   if(my<0){mxI<-(-mx); myI<-(-my)}else{mxI<-mx; myI<-my}
#   
#   if (var.equal == TRUE) {
#     degf <- nx + ny - 2
#     spool <- sqrt((vx * (nx - 1) + vy * (ny - 1))/degf)
#     statistic <- (mxI - myI * rho)/(spool * sqrt(1/nx + (rho^2)/ny))
#     if (alternative == "less") {
#       p.value <- pt(q = statistic, df = degf, lower.tail = TRUE)
#       alpha <- (1 - conf.level)
#     }
#     if (alternative == "greater") {
#       p.value <- pt(q = statistic, df = degf, lower.tail = FALSE)
#       alpha <- (1 - conf.level)
#     }
#     if (alternative == "two.sided") {
#       p.value <- min(1, 2 * pt(q = abs(statistic), df = degf, 
#                                lower.tail = FALSE))
#       alpha <- (1 - conf.level)/2
#     }
#     method <- "Ratio-t-test for equal variances"
#     vpool <- (vx * (nx - 1) + vy * (ny - 1))/degf
#     quant <- qt(p = 1 - alpha, df = degf, lower.tail = TRUE)
#     tA <- ((vpool * quant^2)/ny) - my^2
#     tB <- 2 * mxI * myI
#     tC <- ((vpool * quant^2)/nx) - mx^2
#     if (tA >= 0) {
#       warning("Confidence set unbounded.")
#       upper <- NA
#       lower <- NA
#     }
#     else {
#       upper <- tB/(-2 * tA) - sqrt(((tB/2)^2) - tA * tC)/tA
#       lower <- tB/(-2 * tA) + sqrt(((tB/2)^2) - tA * tC)/tA
#     }
#   }
#   
#   if (var.equal == FALSE & iterativeCI == FALSE) {
#     
#     degf <- max(1, ((vx/nx + (rho^2) * vy/ny)^2)/((vx^2)/((nx^2) * 
#                                                             (nx - 1)) + (rho^4) * (vy^2)/((ny^2) * (ny - 1))) )
#     
#     stderr <- sqrt(vx/nx + (rho^2) * vy/ny)
#     statistic <- (mxI - myI * rho)/stderr
#     
#     if (alternative == "less") {
#       p.value <- pt(q = statistic, df = degf, lower.tail = TRUE)
#       alpha <- (1 - conf.level)
#     }
#     if (alternative == "greater") {
#       p.value <- pt(q = statistic, df = degf, lower.tail = FALSE)
#       alpha <- (1 - conf.level)
#     }
#     if (alternative == "two.sided") {
#       p.value <- min(1, 2 * pt(q = abs(statistic), df = degf, 
#                                lower.tail = FALSE))
#       alpha <- (1 - conf.level)/2
#     }
#     method <- "Ratio t-test for unequal variances"
#     
#     degfest <- max(1, ((vx/nx + (est^2) * vy/ny)^2)/((vx^2)/((nx^2) * 
#                                                                (nx - 1)) + (est^4) * (vy^2)/((ny^2) * (ny - 1))) )
#     
#     quant <- qt(p = 1 - alpha, df = degfest, lower.tail = TRUE)
#     
#     tA <- ((vy * quant^2)/ny) - my^2
#     tB <- 2 * mxI * myI
#     tC <- ((vx * quant^2)/nx) - mx^2
#     
#     if (tA >= 0) {
#       warning("Confidence set unbounded.")
#       upper <- NA
#       lower <- NA
#     }
#     else {
#       upper <- tB/(-2 * tA) - sqrt(((tB/2)^2) - tA * tC)/tA
#       lower <- tB/(-2 * tA) + sqrt(((tB/2)^2) - tA * tC)/tA
#     }
#   }
#   
#   
#   if (var.equal == FALSE & iterativeCI == TRUE) {
#     
#     degf <- ((vx/nx + (rho^2) * vy/ny)^2)/((vx^2)/((nx^2) * 
#                                                      (nx - 1)) + (rho^4) * (vy^2)/((ny^2) * (ny - 1)))
#     stderr <- sqrt(vx/nx + (rho^2) * vy/ny)
#     statistic <- (mxI - myI * rho)/stderr
#     if (alternative == "less") {
#       p.value <- pt(q = statistic, df = degf, lower.tail = TRUE)
#       alpha <- (1 - conf.level)
#     }
#     if (alternative == "greater") {
#       p.value <- pt(q = statistic, df = degf, lower.tail = FALSE)
#       alpha <- (1 - conf.level)
#     }
#     if (alternative == "two.sided") {
#       p.value <- min(1, 2 * pt(q = abs(statistic), df = degf, 
#                                lower.tail = FALSE))
#       alpha <- (1 - conf.level)/2
#     }
#     
#     method <- "Ratio t-test for unequal variances"
#     
#     conf.int <- CIratioiter(nx=nx, ny=ny, mx=mxI, my=myI, vx=vx, vy=vy, alternative = alternative, conf.level = conf.level, ul=ul, ll=ll) 
#     lower<-conf.int[1]
#     upper<-conf.int[2]
#   }
#   
#   
#   if (alternative == "two.sided") {
#     conf.int <- c(lower, upper)
#   }
#   else {
#     if (alternative == "less") {
#       conf.int <- c(-Inf, upper)
#     }
#     else {
#       if (alternative == "greater") {
#         conf.int <- c(lower, Inf)
#       }
#     }
#   }
#   
#   
#   names(statistic) <- "t"
#   estimate <- c(mx, my, est)
#   names(estimate) <- c(paste("mean", namex), paste("mean", 
#                                                    namey), paste(namex, namey, sep = "/"))
#   names(degf) <- "df"
#   names(rho) <- "ratio of means"
#   data.name <- paste(namex, namey, sep = " and ")
#   conf.int<-as.numeric(conf.int)
#   attr(conf.int, "conf.level") <- conf.level
#   
#   out <- list(statistic = statistic, parameter = degf, p.value = p.value, 
#               conf.int = conf.int, estimate = estimate, null.value = rho, 
#               alternative = alternative, method = method, data.name = data.name)
#   class(out) <- "htest"
#   return(out)
# }
# 
