
# MDM Package: most difference in means
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.

# Load package manager
# if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# # p_load(docstring)
# p_load(cubature)


mdm_credint <- function(x, y, conf.level = 0.95, num_param_sims = 250/(1-conf.level), 
                        plot=FALSE, relative = FALSE, sharedVar=FALSE, rand.seed = NA){
  #' @description calculates the (raw/relative) most difference in means, a 
  #' statistic that estimates the largest absolute difference in means supported
  #' by the data. Uses credibility interval.
  #' Citation: https://arxiv.org/abs/2201.01239
  #' 
  #' Equation:
  #' mdm <- Qraw(1-a_dm)
  #' - where Qraw() is the quantile function of the posterior summarizing |u_dm|, 
  #'   and a_dm is the significance level
  #'   u_dm = u_x - u_y
  #'   
  #' rmdm <- Qrel(1-a_dm)
  #' - where Qrel() is the quantile function of the posterior summarizing |ru_dm|, 
  #'   and a_dm is the significance level
  #' 
  #' @param x vector of measurements from group 1 (experiment)
  #' @param y vector of measurements from group 2  (control)
  #' @param conf.level confidence level
  #' @param num_param_sims number of monte carlo trials to calcualte ldm
  #' @param plot plot data
  #' @param relative FALSE: calculates mdm, TRUE: calculates rmdm
  #' @return value of mdm or rmdm

  # save(list = ls(all.names = TRUE), file = "temp/mdm_credint.RData",envir = environment())
  # load(file = "temp/mdm_credint.RData")
  if (!is.na(rand.seed)) {set.seed(rand.seed)}
  
  xbar <- mean(x)
  ybar <- mean(y)
  s2x <- var(x)
  s2y <- var(y)
  m <- length(x)
  n <- length(y)
  if(sharedVar){
    shape <- .5*(m + n - 2)
    scale <- .5*((m-1)*s2x + (n-1)*s2y)
    ssSims <- 1/rgamma(num_param_sims, shape = shape, rate = scale)
    mu1Sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ssSims/m))
    mu2Sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ssSims/n))
  }else{ # different variances
    shape1 <- .5*(m-1)
    scale1 <- .5*(m-1)*s2x
    shape2 <- .5*(n-1)
    scale2 <- .5*(n-1)*s2y
    ss1Sims <- 1/rgamma(n = num_param_sims, shape = shape1, rate = scale1)
    ss2Sims <- 1/rgamma(n = num_param_sims, shape = shape2, rate = scale2)
    mu1Sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ss1Sims/m))
    mu2Sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ss2Sims/n))
  }
  if(!relative){
    cdf <- ecdf(mu1Sims - mu2Sims)
  }else{
    cdf <- ecdf((mu1Sims - mu2Sims)/mu2Sims)
  }
  # solve $F(c) - F(-c) = .95
  upper <- uniroot(function(x){ cdf(x) - cdf(-x) - conf.level},
                   lower = 0,
                   upper = max(c(abs(x),abs(y))),
                   extendInt = "yes")$root
  if(plot & !relative){
    hist(mu1Sims - mu2Sims)
    abline(v=upper,col="red")
    abline(v=-upper,col="red")
  }else if(plot & relative){
    hist((mu1Sims - mu2Sims)/mu2Sims)
    abline(v=upper,col="red")
    abline(v=-upper,col="red")
  }
  return(abs(upper))
}





ldm_credint <- function(x, y, conf.level = 0.95, num_param_sims = 250/(1-conf.level), 
                        plot=FALSE, relative = FALSE, sharedVar=FALSE, keepSign = TRUE, rand.seed = NA){
  #' @description calculates the (raw/relative) least difference in means, a 
  #' statistic that estimates the smallest difference in means supported by the 
  #' data. Uses credibility interval.
  #' 
  #' Equation:
  #' ldm <- sign(x_bar - ybar) * (sign(b_lo) == sign(b_hi)) * max(abs( c(b_lo, b_hi) ))
  #'         original effect sign * force zero if interval includes zero * select closest bound to 0
  #' 
  #' Where b_lo and b_hi are credible bounds for u_dm whjen relative= FALSE and
  #'   r_u_dm when realtive=TRUE
  #' 
  #' @param x vector of measurements from group a, experiment 1
  #' @param y vector of measurements from group a, experiment 1
  #' @param conf.level vector of measurements from group a, experiment 1
  #' @param num_param_sims vector of measurements from group a, experiment 1
  #' @param plot string of number label used for basename of all figures
  #' @param relative path to export figures to disk
  #' @param keepSign base name for exported figures
  #' @return value of ldm ro rldm
  
  # save(list = ls(all.names = TRUE), file = "temp/mdm_credint.RData",envir = environment())
  # load(file = "temp/mdm_credint.RData")
  if (!is.na(rand.seed)) {set.seed(rand.seed)}
  
  xbar <- mean(x)
  ybar <- mean(y)
  s2x <- var(x)
  s2y <- var(y)
  m <- length(x)
  n <- length(y)
  if(sharedVar){
    shape <- .5*(m + n - 2)
    scale <- .5*((m-1)*s2x + (n-1)*s2y)
    ssSims <- 1/rgamma(num_param_sims, shape = shape, rate = scale)
    mu1Sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ssSims/m))
    mu2Sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ssSims/n))
  }else{ # different variances
    shape1 <- .5*(m-1)
    scale1 <- .5*(m-1)*s2x
    shape2 <- .5*(n-1)
    scale2 <- .5*(n-1)*s2y
    ss1Sims <- 1/rgamma(n = num_param_sims, shape = shape1, rate = scale1)
    ss2Sims <- 1/rgamma(n = num_param_sims, shape = shape2, rate = scale2)
    mu1Sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ss1Sims/m))
    mu2Sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ss2Sims/n))
  }
  if(!relative){
    cdf <- ecdf(mu1Sims - mu2Sims)
  }else{
    cdf <- ecdf((mu1Sims - mu2Sims)/mu2Sims)
  }
  
  # TODO use prctile function
  b_lo <- uniroot(function(x){ cdf(x) - conf.level},
                   lower = 0,
                   upper = max(c(abs(x),abs(y))),
                   extendInt = "yes")$root
  
  b_hi <- uniroot(function(x){ cdf(x) - (1-conf.level)},
                  lower = 0,
                  upper = max(c(abs(x),abs(y))),
                  extendInt = "yes")$root

  
  ldm <- sign(xbar - ybar) * (sign(b_lo) == sign(b_hi)) * min(abs( c(b_lo, b_hi) ))
  
  # Keep sign of effect size if requested, since sign matters for practical sig.
  if (!keepSign) {ldm <- abs(ldm)}
  
  # if(plot & !relative){
  #   hist(mu1Sims - mu2Sims)
  #   abline(v=upper,col="red")
  #   abline(v=-upper,col="red")
  # }else if(plot & relative){
  #   hist((mu1Sims - mu2Sims)/mu2Sims)
  #   abline(v=upper,col="red")
  #   abline(v=-upper,col="red")
  # }
  # browser();
  return(ldm)
}





credint <- function(x,y, conf.level= 0.95, num_param_sims = 250/(1-conf.level), 
                    sharedVar=FALSE, relative = FALSE, rand.seed = NA) {
  # x control group
  # y experiment group
  if (!is.na(rand.seed)) {set.seed(rand.seed)}
  
  xbar <- mean(x)
  ybar <- mean(y)
  s2x <- var(x)
  s2y <- var(y)
  m <- length(x)
  n <- length(y)
  if(sharedVar){
    shape <- .5*(m + n - 2)
    scale <- .5*((m-1)*s2x + (n-1)*s2y)
    ssSims <- 1/rgamma(num_param_sims, shape = shape, rate = scale)
    mux_sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ssSims/m))
    muy_sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ssSims/n))
  }else{ # different variances
    shape1 <- .5*(m-1)
    scale1 <- .5*(m-1)*s2x
    shape2 <- .5*(n-1)
    scale2 <- .5*(n-1)*s2y
    ss1Sims <- 1/rgamma(n = num_param_sims, shape = shape1, rate = scale1)
    ss2Sims <- 1/rgamma(n = num_param_sims, shape = shape2, rate = scale2)
    mux_sims  <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ss1Sims/m))
    muy_sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ss2Sims/n))
  }
  if(!relative){
    cdf <- ecdf(muy_sims - mux_sims)
  }else{
    cdf <- ecdf((muy_sims - mux_sims)/mux_sims)
  }
  

  b_lo <- uniroot(function(x){ cdf(x) - (1-conf.level)/2},
                  lower = 0,
                  upper = max(c(abs(x),abs(y))),
                  extendInt = "yes")$root
  
  b_hi <- uniroot(function(x){ cdf(x) - (1-(1-conf.level)/2)},
                  lower = 0,
                  upper = max(c(abs(x),abs(y))),
                  extendInt = "yes")$root
  
  return(c(b_lo, b_hi))
  
}




mdm_confint <- function(x, y, conf.level = 0.95, num_param_sims = 250/(1-conf.level), 
                        plot=FALSE, relative = FALSE){
  # save(list = ls(all.names = TRUE), file = "temp/mdm_credint.RData",envir = environment())
  # load(file = "temp/mdm_credint.RData")
  
  xbar <- mean(x)
  ybar <- mean(y)
  s2x <- var(x)
  s2y <- var(y)
  m <- length(x)
  n <- length(y)
  
  
  mu1Sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(s2x / m))
  mu2Sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(s2y / n))
  # shape1 <- .5*(m-1)
  # scale1 <- .5*(m-1)*s2x
  # shape2 <- .5*(n-1)
  # scale2 <- .5*(n-1)*s2y
  # ss1Sims <- 1/rgamma(n = num_param_sims, shape = shape1, rate = scale1)
  # ss2Sims <- 1/rgamma(n = num_param_sims, shape = shape2, rate = scale2)
  # mu1Sims <- rnorm(n = num_param_sims, mean = xbar, sd = sqrt(ss1Sims/m))
  # mu2Sims <- rnorm(n = num_param_sims, mean = ybar, sd = sqrt(ss2Sims/n))

  # mu1Sims <-rowMeans(matrix(rnorm(n = m*num_param_sims, mean = xbar, sd = sqrt(s2x)), ncol = m))
  # mu2Sims <-rowMeans(matrix(rnorm(n = n*num_param_sims, mean = ybar, sd = sqrt(s2y)), ncol = n))
    
  if(!relative){
    cdf <- ecdf(mu1Sims - mu2Sims)
  } else {
    cdf <- ecdf((mu1Sims - mu2Sims)/mu2Sims)
  }
  # solve $F(c) - F(-c) = .95
  upper <- uniroot(function(x){ cdf(x) - cdf(-x) - conf.level},
                   lower = 0,
                   upper = max(c(abs(x),abs(y))),
                   extendInt = "yes")$root
  if(plot & !relative){
    hist(mu1Sims - mu2Sims)
    abline(v=upper,col="red")
    abline(v=-upper,col="red")
  }else if(plot & relative){
    hist((mu1Sims - mu2Sims)/mu2Sims)
    abline(v=upper,col="red")
    abline(v=-upper,col="red")
  }
  return(abs(upper))
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
  
  # Quantile with prob set to alpha and sample estimates as parameters
  mdm <-  tryCatch(
    qft(p = conf.level, df = df_d, mu = xbar_dm, sigma = sd_dm,
        lo.bound, hi.bound),
    error = function(c) NaN)
  
  # If integration fails pause execution
  if (is.nan(mdm) || is.null(mdm)) {browser();}
  
  return(mdm)
}



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
  df <- n_ctrl + n_exp - 2
  
  sSquared <- (sum((x_ctrl - mean_ctrl)^2) + sum((y_exp - mean_exp)^2))/df
  
  fieller_int <- tryCatch(get_FiellerInterval(mean_ctrl, mean_exp, sSquared, 
                                              1/n_ctrl, 1/n_exp, df, v12 = 0, alpha=1-conf.level), 
                          error = function(c) data.frame(upper = NaN, lower = NaN))
  
  rmdm <- max(abs(c(fieller_int$lower,fieller_int$upper)))
  
  
  return(rmdm)
}





ldm_tdist <- function(x, y = NULL, conf.level = 0.95) {
  #' @description Calculate least difference in means statistic from integrating 
  #' over a folded t-distribution to an extent dictated by (1-conf.level)
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @return Returns most difference in means (single value)
  #' @usage mdm_tdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=6,mean=0,sd=1); y <- rnorm(n=6,mean=0,sd=1);
  #' ldm_tdist(x,y, conf.level = 0.95)
  
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
    wtt <- t.test(x=x,y=y, var.equal = FALSE, conf.level = 1-(1-conf.level)/4)
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
  
  # Calculate search bounds for ldm, used for uniroot call
  lo.bound <- 0
  hi.bound <- abs(xbar_dm)
  if (lo.bound>hi.bound) {browser();}
  
  
  # Calculate lower bound of mu_dm
  lo.bound <- abs(xbar_dm) - qt(conf.level, df=df_d)*sd_dm
  # Equivalent to:
  # wtt <-t.test(x=y,y=x, var.equal = FALSE, conf.level = 1-2*(1-conf.level))
  if (lo.bound < 0) {
    ldm <- 0
  } else {
    # Lower quantile of folded normal
    ldm <-  tryCatch(qft(p = 1-conf.level, df = df_d, mu = xbar_dm, sigma = sd_dm,
                         lo.bound, hi.bound), error = function(c) NaN)
  }
  
  # If integration fails pause execution
  if (is.nan(ldm) || is.null(ldm)) {browser();}
  
  return(ldm)
}



lacb_tdist_2sample <- function (x, y, conf.level = 0.95) {
  #' Least absolute one-tailed confidence bounds of t-distribution
  #' 
  lo_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "greater")$conf.int[1]
  up_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "less")$conf.int[2]
  
  lacb <- min(abs(c(lo_b, up_b)))
  if (sign(lo_b)!= sign(up_b)) {lacb <- 0}
  
  return(lacb)
}


rldm_tdist <- function(x_ctrl, y_exp, conf.level = 0.95, 
                       verbose = FALSE,  var.equal = FALSE)  {
  #' @description Calculates the relative most difference in means assuming in 
  #' the form of rldm = (y_exp - x_ctrl)/x_ctrl
  #' 
  #' @param x_ctrl vector of measurements in control group
  #' @param y_exp vector of measurements in experimental group
  #' @param conf.level significance level for calculating ldm
  #' 
  #' @return relative least difference in means
  # Equation for pooled variance taken from:
  # https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals5.html
  
  # Calculate basic sample statistics
  mean_ctrl = mean(x_ctrl); sd_ctrl = sd(x_ctrl)
  n_ctrl = length(x_ctrl);se_ctrl = sd_ctrl/sqrt(n_ctrl)
  
  mean_exp = mean(y_exp); sd_exp = sd(y_exp)
  n_exp = length(y_exp); se_exp = sd_exp/sqrt(n_exp)
  
  mean_dm = mean_exp - mean_ctrl
  df <- n_ctrl + n_exp - 2
  
  sSquared <- (sum((x_ctrl - mean_ctrl)^2) + sum((y_exp - mean_exp)^2))/df
  
  fieller_int <- tryCatch(get_FiellerInterval(mean_ctrl, mean_exp, sSquared, 
                                              1/n_ctrl, 1/n_exp, df, v12 = 0, alpha=1-conf.level), 
                          error = function(c) data.frame(upper = NaN, lower = NaN))
  
  # RLDM is the bounds closest to zero, or zero of the bounds flank zero
  rldm <- min(abs(c(fieller_int$lower,fieller_int$upper))) * 
    (sign(fieller_int$lower)== sign(fieller_int$upper))
  
  return(rldm)
}




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
  # save(list = ls(all.names = TRUE), file = "temp/qft.RData",envir = environment())
  # load(file = "temp/qft.RData")
  
  # Integration is also used to calculate p, but it can be a very sparse 
  # integration (long spans where f(x)=0 and then a very small width spike in the pdf)
  # To deal with this, integration is broken up into separate intervals to make 
  # sure the spike is not missed
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
  
  return(xroot)
}

get_FiellerInterval <- function(xbar, ybar, sSquared, v11, v22, f, v12 = 0, alpha=.025){
  tQuantile <- qt(1-alpha, f)
  A <- xbar^2 - v11*sSquared*tQuantile^2
  if(A <= 0)
    stop("confidence interval not available (unless you are okay with two disjoint intervals)")
  B <- -2*((ybar - xbar)*xbar + sSquared*tQuantile^2*(v11 - v12))
  C <- (ybar - xbar)^2 - sSquared*tQuantile^2*(v22 + v11 - 2*v12)
  discriminant <- B^2 - 4*A*C
  if(discriminant <= 0)
    stop("confidence interval not available (complex-valued)")
  center <- -B/2/A
  width <- sqrt(discriminant)/2/A
  list(lower = center - width, upper = center + width)
}

