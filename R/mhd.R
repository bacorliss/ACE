#(broom)


# MHD: Most Hidden Difference
# Calculates the upper bounds of effecct size from a single distribution
# or difference between two distributions.




t_dist_2s_unpaired <- function(x, y, alpha = 0.05) {
  nx = length(x)
  ny = length(y)
  d_bar_ppoled = mean(x)-mean(y)
  df_pooled = nx+ny-2
  sd_pooled = ((nx-1)*sx^2 + (ny-1)*sy^2)/(nx+ny-2)
  
  ucl = t_dist(x_bar = d_bar_ppoled, sx = sd_pooled, nx = nx, dfx = dfx, alpha = 0.05)
  
  return(ucl)
}

t_dist_2s_paired <- function(x, y, alpha = 0.05) {
  # Pairwise difference between observations
  d = x-y
  ucl = t_dist_1s(d, alpha)
  return(ucl)
}

t_dist_1s <- function(x, alpha = 0.05) {
# Calculates upper bounds of mean of folded t-distribution based on pre-folded 
# central t-distribution
  
  # Sample statistics
  x_bar <- mean(x)
  sx <- sd(x)
  nx <- length(x)
  dfx <- n - 1
  
  ucl = folded_t_dist(x_bar = x_bar, sx = sx, nx = nx, dfx = dfx, alpha = 0.05)
  
  return(ucl)
}


folded_t_dist <- function(x_bar, sx, nx, dfx, alpha = 0.05) {
  # Estimate multiplier for bounds of search interval
  sd_mult <- qt(1 - alpha, dfx)
  
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  t_star_function <- function(x) {abs(pt(-x_bar + x,dfx) - pt(-x_bar - x, dfx) - (1 - alpha))}
  
  # Calculate roughly the extext of search space
  search_bounds = c(0, x_bar + 4*sd_mult * s)
  
  # Minimization
  t_star = optimize(t_star_function, search_bounds, maximum = FALSE)
  
  if (t_star$minimum == search_bounds[2]) {
    warning("mhd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # CI_mean - x_bar +- t_star * s/sqrt(n)
  ucl = x_bar + t_star$minimum * s/sqrt(n)
  
  return(ucl)
}


z_dist_1s <- function(x, alpha = 0.05) {
  # Calculates upper bounds of mean of folded t-distribution based on pre-folded 
  # central t-distribution
  
  # Sample statistics
  x_bar <- mean(x)
  s <- sd(x)
  n <- length(x)
  df <- n - 1
  
  # Calculate upper confidence limit of mean of folded distribition
  ucl = folded_z_dist(x_bar, sx, nx, dfx, alpha)
    
    return(ucl)
}
 

folded_z_dist <- function(x_bar, sx, nx, dfx, alpha = 0.05) { 
  # Estimate multiplier for bounds of search interval
  sd_mult <- qt(1 - alpha, dfx)
  
  # Integration of folded t-distribution can be calculate from standard central t-distribution
  z_star_function <- function(x) {abs(pnorm(-x_bar + x,sx) - pnorm(-x_bar - x, sx) - (1 - alpha))}
  
  # Calculate roughly the extext of search space
  search_bounds = c(0, x_bar + 4*sd_mult * sx)
  
  # Minimization
  z_star = optimize(z_star_function, search_bounds, maximum = FALSE)
  
  if (z_star$minimum == search_bounds[2]) {
    warning("mhd: endpoint minimization equal to upper bounds of search space. Results unreliable.")
  }
  
  # CI_mean - x_bar +- z_star * s/sqrt(n)
  ucl = x_bar + z_star$minimum * x/sqrt(n)

  return(ucl)  
}
