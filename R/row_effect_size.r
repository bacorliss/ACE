


# USer defined functions
# Calculate variance by row like rowMeans or rowSums
# sum(x-x_bar).^2/(n-1)
rowVars <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)}
rowSS <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)}

compute_effect_size <- function(x,y) {
  
  # Stats of x
  means_x = rowMeans(x)
  ss_x = rowSS(x)
  n_x = dim(x)[2]
  vars_x = ss_x / n_x
  
  # Stats of y
  means_y = rowMeans(y)
  ss_y = rowSS(y)
  n_y = dim(y)[2]
  vars_y = ss_y / n_y
  
  # Pooled standard deviation
  s_pooled = sqrt(  ( (n_x - 1) * vars_x + (n_y - 1) * vars_y) /
                    (n_x + n_y - 2) )
  
  # Approx. correction factor for hedges g
  g_star = (1 - 3/ (4 * (n_x + n_y) - 9)  )
  
  # Difference between means of two distributuons
  df$difference_of_means = x_mean - y_mean
  # Difference in variance between two distributions
  df$difference_of_variances = x_var-y_var
  # Difference in means normalized to standard deviation
  df$cohens_d = (means_x-means_y)/s_pooled
  # Difference in means normalized to the standard deviation of the control sample
  df$glass_delta = (means_x-means_y) / sqrt(vars_y)
  # Differences in means normalized to pooled standard deviation and then corrected
  df$hedges_g = g_star * (means_x-means_y) / s_pooled

}