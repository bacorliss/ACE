
# MMD Package: Most Mean Difference
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.



mmd_normal(x, y = NULL, paired = FALSE, var.equal = FALSE, 
           conf.level = 0.95, verbose = FALSE, distribution = NULL){
  if (is.null(x)) errorCondition('input argument x must be specified.')
  
  # Determine if data is 1-sample of 2-sample, calcaulte mmd
  if (is.null(y) || paired) {
    # Convert 2 sample paired to 1 sample
    if (paired) d=y-x else d=x
    mmd <-  mmd_1_sample_normal(x = d, var.equal = var.equal, 
                               conf.level = conf.level, verbose = verbose, 
                               distribution = distribution)
  } else {
    mmd <-  mmd_2_sample_normal(x = x, y = y, paired = paired, var.equal = var.equal, 
                                conf.level = conf.level, verbose = verbose, 
                                distribution = distribution)
  }
  return(mmd)
}


mmd_1sample_normal <- function(x, conf.level, distribution = NULL) {
  #  Calculates mmd for 1 sample, defined as the upper confidence limit of mean for 1 sample 
  # based on a folded t-distibution
  # Input
  #   x: sample of observations
  #   alpha: confidence level
  #   distribution: distribution distribution, either t-distribution or z-distribution, with NULL
  #              it is chosen based on whether df > 30.
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Sample statistics
  alpha <- 1-conf.level
  x_bar <- mean(x)
  sx <- sd(x)
  nx <- length(x)
  dfx <- nx - 1
  
  # Determine distribution if not specified
  if (dfx<30-1 & is.null(distribution)) {distribution <- 't_dist'} else {distribution<-'z_dist'}
  distribution <- 't_dist'
  
  if (distribution == 't_dist'){
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_tdist(x_bar = x_bar, sx = sx, nx = nx, dfx = dfx, alpha = 0.05)
  } else if (distribution == 'z_dist') {
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_zdist(x_bar = x_bar, sx = sx, nx = nx, dfx = dfx, alpha = 0.05)
  }
  return(ucl)
}


mmd_2sample_normal <- function(x, y, alpha = 0.05, distribution = NULL) {
  #  Calculates mmd for two unpaired samples, defined as the upper confidence limit of the mean
  # based on a folded t-distibution
  # Input
  #   x: sample of observations
  #   y: sample of observations
  #   alpha: confidence level
  #   distribution: distribution distribution, either t-distribution or z-distribution, with NULL 
  #              it is chosen based on whether df > 30.
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Sample number and standard deviations
  nx <- length(x); ny <- length(y)
  sx <- sd(x); sy <-  sd(y)
  # Pooled mean, degrees of freedom, and standard deviation
  d_bar <- mean(x) - mean(y)
  df_d_pooled <- nx + ny - 2
  sd_d_pooled <- ((nx-1)*sx^2 + (ny-1)*sy^2)/df
  
  #  Standard error of the mean of d: sd multipled by sqrt of the sum of 
  # reciprocal squares of n for each group
  sem_pooled <- sd_pooled * sqrt((1 / nx^2) + (1 / ny^2))
  
  # Determine distribution if not specified
  if (dfd_pooled<30-2 & is.null(distribution)) {distribution <-'t_dist'} else {distribution <- 'z_dist'}
  distribution <- 't_dist' #TODO restore IF statement
  
  if (distribution == 't_dist'){
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_tdist(x,y, alpha = 0.05)
  } else if (distribution == 'z_dist') {
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_zdist(x_bar = d_bar, sx = sd_pooled, nx = nx, #TODO This is wrong
                                semx = semd_pooled, dfx = df_pooled, alpha = 0.05)
  }
  return(ucl)
}




mmd_2sample_normal_tdist <- function(x,y, conf.level = 0.95, verbose = FALSE, search_pad_percent = 0.1) {
  #  Calculate t statistic from integrating a central normal pdf shifted by -x_bar
  # Using minimization function to integrate over a normal CDF with area 
  # under the curve to (1-a)
  # Inputs
  #   x_bar: sample mean
  #   sx: Sample standard deviation
  #   nx: sample number of first grouo
  #   
  #   dfx: sample degrees of freedom
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Sample number and standard deviations
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  # Pooled mean, degrees of freedom, and standard deviation
  d_bar <- mean(x) - mean(y)
  df_d <- n_x + n_y - 2
  sd_d <- ( (n_x-1)*sd_x^2 + (ny - 1) * sd_y^2 ) / df
  sem_d = sd_d / sqrt(nx) # TODO wrong
  df_d = nx + ny - 2
  
  # Calculate search bounds defined by tails of alpha and 2*alpha CI of mean 
  ci_mean_alpha <- t.test(x,y, conf.level = conf.level, paired = FALSE, 
                       var.equal = FALSE, alternative = "two.sided")$conf.int
  ci_mean_2alpha <- t.test(x,y, conf.level = 1-2*(1-conf.level), paired = FALSE, 
                       var.equal = FALSE, alternative = "two.sided")$conf.int
  lower_bounds =  min(abs(ci_mean_alpha))
  upper_bounds = max(abs(ci_mean_2alpha))
  bounds_range = upper_bounds - lower_bounds
  search_bounds = c(min(c(lower_bounds - search_pad_percent * bounds_range, 0)),
                    upper_bounds + search_pad_percent * bounds_range)
  
  # Adjust the obs. # for t score normalzation (2 sample or 1 sample)
  if (!is.null(ny))  t_score_norm <- sqrt(1/nx+1/ny) else t_score_norm <- sqrt(1/nx)
  
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  t_star_function <- function(x) {pt(-d_bar / (sd_d * t_score_norm) + x, df_d) - 
                                  pt(-d_bar / (sd_d * t_score_norm) - x, df_d) - conf.level}
  
  # Solve for t star with root finding where t_star_function equals (1-alpha)
  t_star = uniroot(t_star_function, search_bounds, check.conv = TRUE,
                         tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
  
  # The optimized root should fall entirely within the earch bounds 
  search_bounds_check <- function(t_star, search_bounds, verbose = FALSE, range_tol=1000)
    
  # CI_mean - x_bar +- t_star * s/sqrt(n)
  ucl = abs(x_bar) + (t_star$root - abs(x_bar) / (sd_d*t_score_norm)) * sem_d
  
  return(ucl)
}


search_bounds_check <- function(t_star, search_bounds, verbose = FALSE, range_tol=1000) {
  
  # Plot root finding data and root if verbose mode on
  if (verbose) {
    x <- seq(search_bounds[1], search_bounds[2], diff(search_bounds)/100)
    y <- lapply(x,t_star_function)
    plot(x,y)
    print(sprintf("x_bar: %.2f, t_star: %.2f",x_bar, t_star$root))
  }
  
  if (abs(t_star$root - search_bounds[2]) < abs(diff(search_bounds)/range_tol) ) {
    warning("mmd: mmd equation root equal to upper bounds of search space: 
            results unreliable.")
  }
  
  if (abs(t_star$root - search_bounds[1]) < abs(diff(search_bounds)/range_tol)  ) {
    warning("mmd: mmd equation root equal to lower bounds of search space: 
            results unreliable.")
    }
  
  
}

ucl_mean_folded_zdist <- function(x_bar, sx, nx, dfx, semx = NULL, alpha = 0.05) { 
  #  Calculate z statistic from integrating a central t-distribution pdf shifted by -x_bar
  # Using minimization function to integrate using a t-distribution CDF setting the area 
  # under the curve to (1-a)
  # Inputs
  #   x_bar: sample mean
  #   sx: Sample standard deviation
  #   nx: sample number
  #   dfx: sample degrees of freedom
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Estimate multiplier for bounds of search interval
  sd_mult <- qnorm(1 - alpha, dfx)
  
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  z_star_function <- function(x) {abs(pnorm(-x_bar/sx + x,sx) - pnorm(-x_bar/sx - x, sx) - (1 - alpha))}
  
  # Calculate roughly the extext of search space
  search_bounds = c(abs(x_bar)/sx, abs(x_bar)/sx + 5*sd_mult)
  
  # Solve for t star with root finding
  t_star = uniroot(t_star_function, search_bounds, check.conv = TRUE,
                   tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
  
  if (z_star$minimum == search_bounds[2]) {
    warning("mmd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = sx / sqrt(nx)}
  
  # CI_mean - x_bar +- z_star * s/sqrt(n)
  ucl = z_star$minimum * x/sqrt(x)

  return(ucl)  
}








ucl_tdist_mean <- function(x_bar, sx, nx, dfx, semx = NULL, alpha = 0.05) {
  #  Calculate t statistic from integrating a central normal pdf shifted by -x_bar
  # Using minimization function to integrate using a normal CDF setting the area 
  # under the curve to (1-a)
  # Inputs
  #   x_bar: sample mean
  #   sx: Sample standard deviation
  #   nx: sample number
  #   dfx: sample degrees of freedom
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Estimate multiplier for bounds of search interval
  sd_mult <- qt(1 - alpha, dfx)
  
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  t_star_function <- function(x) { abs( (pt(x,dfx) - pt(-x, dfx))  - (1 - alpha))}
  
  # Calculate roughly the extext of search space
  search_bounds = c(0, x_bar + 4*sd_mult * sx)
  
  # Minimization
  t_star = optimize(t_star_function, search_bounds, maximum = FALSE)
  
  if (t_star$minimum == search_bounds[2]) {
    warning("mmd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = sx / sqrt(nx)}
  
  # CI_mean x_bar +- t_star * s/sqrt(n)
  ucl = x_bar + t_star$minimum * semx
  
  print(t_star$minimum)
  
  return(ucl)
}

