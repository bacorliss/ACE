
# MHD: pakage: Most Hidden Difference
#  Calculates the upper bounds of the mean of the effect size from a single distribution
# or difference between two distributions.


ucl_mean_tdist <- function(x_bar, sx, nx, dfx, semx = NULL, alpha = 0.05) {
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
    warning("mhd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = sx / sqrt(nx)}
  
  # CI_mean x_bar +- t_star * s/sqrt(n)
  ucl = x_bar + t_star$minimum * semx
  
  print(t_star$minimum)
  
  return(ucl)
}


mhd_1sample <- function(x, alpha = 0.05, dist_type = NULL) {
  #  Calculates MHD for 1 sample, defined as the upper confidence limit of mean for 1 sample 
  # based on a folded t-distibution
  # Input
  #   x: sample of observations
  #   alpha: confidence level
  #   dist_type: distribution dist_type, either t-distribution or z-distribution, with NULL
  #              it is chosen based on whether df > 30.
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  # Sample statistics
  x_bar <- mean(x)
  sx <- sd(x)
  nx <- length(x)
  dfx <- nx - 1
  
  # Determine distribution if not specified
  if (dfx<30-1 & is.null(dist_type)) {dist_type <- 't_dist'} else {dist_type<-'z_dist'}

  if (dist_type=='t_dist'){
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_tdist(x_bar = x_bar, sx = sx, nx = nx, dfx = dfx, alpha = 0.05)
  } else if (dist_type=='z_dist') {
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_zdist(x_bar = x_bar, sx = sx, nx = nx, dfx = dfx, alpha = 0.05)
  }

  return(ucl)
}


mhd_2sample_paired <- function(x, y, alpha = 0.05, dist_type = NULL) {
  #  Calculates MHD for two paired samples, defined as the upper confidence limit of the mean
  # based on a folded t-distibution
  # Input
  #   x: sample of observations
  #   y: sample of observations
  #   alpha: confidence level
  #   dist_type: distribution dist_type, either t-distribution or z-distribution, with NULL 
  #              it is chosen based on whether df > 30.
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  # Effect distribution is sample by sample different between x and y
  d = x-y
  
  # Two sample paired is equivalent to 1 sample
  ucl = mhd_1sample(d, alpha, dist_type=dist_type)
  return(ucl)
}


mhd_2sample_unpaired <- function(x, y, alpha = 0.05, dist_type = NULL) {
  #  Calculates MHD for two unpaired samples, defined as the upper confidence limit of the mean
  # based on a folded t-distibution
  # Input
  #   x: sample of observations
  #   y: sample of observations
  #   alpha: confidence level
  #   dist_type: distribution dist_type, either t-distribution or z-distribution, with NULL 
  #              it is chosen based on whether df > 30.
  # Output
  #   ucl: upper confidence limit of the mean of folded distribution
  
  nx <- length(x)
  ny <- length(y)
  # Pooled mean, degrees of freedom, and standard deviation
  d_bar_pooled <- mean(x)-mean(y)
  dfd_pooled <- nx + ny - 2
  sd_pooled <- ((nx-1)*sx^2 + (ny-1)*sy^2)/df
  
  #  Standard error of the mean of d: sd multipled by sqrt of the sum of 
  # reciprocal squares of n for each group
  semd_pooled <- sd_pooled * sqrt((1 / nx^2) + (1 / ny^2))
  
  # Determine distribution if not specified
  if (dfd_pooled<30-2 & is.null(dist_type)) {dist_type<-'t_dist'} else {dist_type<-'z_dist'}
  
  if (dist_type=='t_dist'){
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_tdist(x_bar = d_bar_pooled, sx = sd_pooled, 
                                semx = semd_pooled, dfx = dfd, alpha = 0.05)
  } else if (dist_type=='z_dist') {
    # Calculate upper confidence limit of mean of folded normal
    ucl = ucl_mean_folded_zdist(x_bar = d_bar_pooled, sx = sd_pooled, 
                                semx = semd_pooled, dfx = dfd, alpha = 0.05)
  }
  
  return(ucl)
}


ucl_mean_folded_tdist <- function(x_bar, sx, nx, dfx, semx = NULL, alpha = 0.05) {
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
  t_star_function <- function(x) {abs(pt(-x_bar / (sx/sqrt(nx)) + x, dfx) - 
                                      pt(-x_bar / (sx/sqrt(nx)) - x, dfx) - (1 - alpha))}
  
  # Calculate roughly the extext of search space
  search_bounds = c(abs(x_bar) / (sx/sqrt(nx)), abs(x_bar) / (sx/sqrt(nx)) + 5*sd_mult)
  #print(sprintf('USL:%f',abs(x_bar) / (sx/sqrt(nx)) + 4*sd_mult*sx))
  
  
  x <- seq(search_bounds[1], search_bounds[2], diff(search_bounds)/100)
  y <- lapply(x,t_star_function)
  plot(x,y)
  
  
  # Minimization
  t_star = optimize(t_star_function, search_bounds, maximum = FALSE, tol=.Machine$double.eps^(1/2))
  # t_star2 = uniroot(t_star_function, search_bounds, check.conv = FALSE,
  #         tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
  # print(t_star2)
  # 
  print(sprintf("x_bar: %.2f, t_star: %.2f",x_bar, t_star$minimum))
  
  if (abs(t_star$minimum - search_bounds[2]) < diff(search_bounds)/1000 ) {
    warning("mhd: minimized point equal to max bounds of search space- results unreliable.")
  }

  if (abs(t_star$minimum - search_bounds[1]) < sx/2 ) {
    warning("mhd: minimized point minimization equal to min bounds of search space- results unreliable.")
  }

  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = sx / sqrt(nx)}
  
  # CI_mean - x_bar +- t_star * s/sqrt(n)
  ucl = abs(x_bar) + (t_star$minimum - abs(x_bar)/(sx/sqrt(nx))) * semx
  
  return(ucl)
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
  search_bounds = c(0, x_bar + 4*sd_mult * sx)
  
  # Minimization
  z_star = optimize(z_star_function, search_bounds, maximum = FALSE)
  
  if (z_star$minimum == search_bounds[2]) {
    warning("mhd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # If SEM is not specified, calculate
  if (is.null(semx)) { semx = sx / sqrt(nx)}
  
  # CI_mean - x_bar +- z_star * s/sqrt(n)
  ucl = z_star$minimum * x/sqrt(n)

  return(ucl)  
}
