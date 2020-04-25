
# MMD Package: Most Mean Difference
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.

# TODO var.equal = TRUE not coded yet

library(docstring)



mmd_normal <- function(x, y = NULL, paired = FALSE, var.equal = FALSE, conf.level = 0.95, 
           verbose = FALSE, distribution = NULL) {
  #' Calculate Most Mean Difference using a normal distribution (z or t)
  #'
  #' @description  Most mean difference is the largest 
  #' difference than could exist between group(s) in the data.
  #' 
  #' Calculates most mean difference assuming a normal distribution
  #' (z or t) by integrating a central normal pdf shifted by sample mean and
  #' Using a root finding function. Most mean difference is the largest 
  #' difference than could exist between group(s) in the data.
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param paired boolean that specifies if 2-sample data is paired (default: FALSE)
  #' @param var.equal boolean if variance assumed equal (default: FALSE)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @param verbose function prints out extra input during execution
  #' @param distribution string specifying assumed distribution ('t_dist', 
  #' 'z_dist') (Default: NULL). If not specified, t_dist chosen for cases under 
  #' 30 samples
  #' @return Returns most mean difference (single value)
  #' @usage mmd_normal_tdist(x)
  #' @usage mmd_normal_tdist(x, y)
  #' @usage mmd_normal_tdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=10,mean=0,sd=1); 
  #' y <- rnorm(n=10,mean=0,sd=1); 
  #' mmd_normal(x,y)
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
  
  # Calculate MMD with specified distribution
  if (distribution=='t_dist'){
    mmd <- mmd_normal_tdist(x, y, conf.level = conf.level, verbose = verbose,
                                     search_pad_percent = 0.1,var.equal = var.equal)
  } else if (distribution =='z_dist') {
    mmd = mmd_normal_zdist(x,y, conf.level = conf.level, verbose = verbose,
                                search_pad_percent = 0.1,var.equal = var.equal)
  } else {
    mmd=NULL
    warning("Input argument 'distribution' not supported")
  }
  
  return(mmd)
}


mmd_normal_tdist <- function(x,y = NULL, conf.level = 0.95, verbose = FALSE, 
                             var.equal = FALSE, search_pad_percent = 0.25) {
  #' Calculate Most Mean Difference using t distribution
  #'
  #' @description  Calculate most mean difference statistic from integrating a
  #' central normal t-distribution pdf shifted by -x_bar Using a root finding 
  #' function to integrate over a normal CDF with area  under the curve equal 
  #' to (1-a). Calculate stats of the difference in means distribution for two 
  #' samples, or keeps same distribution for one sample case.
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @param verbose function prints out extra input during execution
  #' @param search_pad_percent for calculating statistic, specifies how far 
  #' outside the theoretical search window the root finding can look.
  #' @return Returns most mean difference (single value)
  #' @usage mmd_normal_tdist(x)
  #' @usage mmd_normal_tdist(x, y)
  #' @usage mmd_normal_tdist(x, y, conf.level = 0.95)
  #' @examples 
  #' x <- rnorm(n=8,mean=0,sd=1); 
  #' y <- rnorm(n=8,mean=0,sd=1); 
  #' mmd_normal_tdist(x,y)

  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distirbution of the sampel for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample stats
    d_bar <- mean(x)
    df_d <- n_x - 1
    sd_d <- sd_x
    sem_d = sqrt( sd_x^2 / n_x)
    
  } else { 
    # 2-sample stats
    d_bar <- mean(y) - mean(x)
    df_d <- n_x + n_y - 2
    sd_d <- sqrt( (( n_x-1)*sd_x^2 + (n_y - 1) * sd_y^2 ) / df_d)
    sem_d = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
    if (verbose) print(sprintf("d_bar: %.3f", d_bar))  
  }
  
  # Calculate search bounds defined by tails of alpha and 2*alpha CI of mean 
  #alpha = (1 - conf.level)
  #ci_mean_alpha  <- qt( c(  alpha/2, 1 -   alpha/2), ncp = d_bar, sd = sd_d)
  #ci_mean_2alpha <- qt( c(2*alpha/2, 1 - 2*alpha/2), ncp = d_bar, sd = sd_d)
  
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
  
  # Calculate MMD with root finding optmization
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  t_star_standard <- function(x) {pt(q = (-d_bar + x) / sem_d, df = df_d) - 
                                  pt(q = (-d_bar - x) / sem_d, df = df_d) - conf.level}
  
  # Debugging search bounds are correct
  if (verbose) {
    t = seq(from = search_bounds[1], to = search_bounds[2], by = diff(search_bounds)/100)
    f_t = sapply(t, t_star_standard)
    plot(t,f_t)
  }
   
  # Solve for t star with root finding where t_star_standard equals (1 - alpha)
  t_star = uniroot(t_star_standard, search_bounds, check.conv = TRUE,
                         tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
  # The optimized root should fall entirely within the earch bounds 
  search_bounds_check(t_star, search_bounds, verbose = verbose, range_tol = 1000)
    
  # MMD is root location added to difference of means
  if (verbose) print(sprintf("t_star: %.3f", t_star$root))
  mmd = t_star$root 
  
  return(mmd)
}

mmd_normal_zdist <- function(x, y = NULL, conf.level = 0.95, verbose = FALSE, 
                             var.equal = FALSE, search_pad_percent = 0.1) {
  #' Calculate Most Mean DIfference using z distribution
  #'
  #' @description Calculate most mean difference statistic from integrating a 
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
  #' @return Returns most mean difference (single value)
  #' @usage mmd_normal_zdist(x)
  #' @usage mmd_normal_zdist(x, y)
  #' @usage mmd_normal_zdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=50,mean=0,sd=1); y <- rnorm(n=50,mean=0,sd=1); 
  #' mmd_normal_zdist(x,y)
  
  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distirbution of the sampel for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample stats
    df_d <- n_x - 1
    d_bar <- mean(x)
    sd_d <- sd_x
    sem_d = sqrt( sd_x^2 / n_x)
    
  } else { 
    # 2-sample stats
    df_d <- n_x + n_y - 2
    d_bar <- mean(y) - mean(x)
    sd_d <- sqrt( (( n_x - 1) * sd_x^2  +  (n_y - 1) * sd_y^2 ) / df_d)
    sem_d = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
    if (verbose) print(sprintf("d_bar: %.3f", d_bar))  
  }
  
  # Calculate search bounds defined by tails of alpha and 2*alpha CI of mean 
  alpha = (1 - conf.level)
  ci_mean_alpha  <- qnorm( c(  alpha/2, 1 -   alpha/2), mean = d_bar, sd = sem_d)
  ci_mean_2alpha <- qnorm( c(2*alpha/2, 1 - 2*alpha/2), mean = d_bar, sd = sem_d)
  lower_bounds =  max(abs(ci_mean_2alpha))
  upper_bounds = max(abs(ci_mean_alpha))
  # Add extra padding around search bounds so root finding not done on boundary
  bounds_range = upper_bounds - lower_bounds
  search_bounds = c(lower_bounds - search_pad_percent * bounds_range,
                    upper_bounds + search_pad_percent * bounds_range)
  if (verbose) print(sprintf('Bounds:[ %.3f  %.3f]', search_bounds[1], search_bounds[2]))
  
  # Integration of folded z-distribution from standard central z-distribution
  z_star_noncentral <- function(x) {pnorm( +x, mean = d_bar, sd = sem_d) - 
    pnorm( -x, mean = d_bar, sd = sem_d) - conf.level}
  # Note: d_bar and x do not need to be in z_score units because they are z 
  # normalized after insertion into function
  z_star_standard <- function(x) { pnorm( (-d_bar + x)/sem_d, mean = 0, sd = 1) - 
    pnorm( (-d_bar - x)/sem_d, mean = 0, sd = 1) - conf.level}
  
  if (verbose) {
    z = seq(from = search_bounds[1], to = search_bounds[2], by = diff(search_bounds)/100)
    f_z = lapply(z, z_star_standard)
    plot(z,f_z)
  }
  
  # Solve for t star with root finding
  z_star = uniroot(z_star_standard, search_bounds, check.conv = TRUE,
                   tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
  # The optimized root should fall entirely within the earch bounds 
  search_bounds_check(z_star, search_bounds, verbose = FALSE, range_tol=1000)
  
  # MMD is root of integration
  if (verbose) print(sprintf("t_star: %.3f", z_star$root))
  mmd = z_star$root
  
  return(mmd)  
}



search_bounds_check <- function(t_star, search_bounds, verbose = FALSE, range_tol=1000) {
  
  
  if (abs(t_star$root - search_bounds[2]) < abs(diff(search_bounds)/range_tol) ) {
    warning("mmd: equation root equal to upper bounds of search space: 
            results unreliable.")
  }
  
  if (abs(t_star$root - search_bounds[1]) < abs(diff(search_bounds)/range_tol)  ) {
    warning("mmd: equation root equal to lower bounds of search space: 
            results unreliable.")
  }
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
    warning("mmd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = s_x / sqrt(n_x)}
  
  # CI_mean x_bar +- t_star * s/sqrt(n)
  ucl = x_bar + t_star$minimum * sem_x
  
  print(t_star$minimum)
  
  return(ucl)
}

