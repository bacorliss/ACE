




lacb_tdist_2sample <- function (x, y, conf.level = 0.95) {
  #' Least absolute one-tailed confidence bounds of t-distribution
  #' 
  lo_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "greater")$conf.int[1]
  up_b <- t.test(x = x, y = y, conf.level = conf.level, alternative = "less")$conf.int[2]
  
  lacb <- min(abs(c(lo_b, up_b)))
  if (sign(ci[1])!=signci[2]) {lacb <- 0}
  
  return(lacb)
}



ldm_normal <- function(x, y = NULL, paired = FALSE, var.equal = FALSE, conf.level = 0.95, 
                       verbose = FALSE, distribution = NULL) {
  #' @description  most difference in means is the largest 
  #' difference than could exist between group(s) in the data.
  #' 
  #' Calculates the least difference in means based on the confidence limits
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
  #' @usage ldm_normal_tdist(x)
  #' @usage ldm_normal_tdist(x, y)
  #' @usage ldm_normal_tdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=10,mean=0,sd=1); 
  #' y <- rnorm(n=10,mean=0,sd=1); 
  #' ldm_normal(x,y)
  #' 
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
  
  # Calculate ldm with specified distribution
  if (distribution=='t_dist'){
    ldm <- ldm_normal_tdist(x, y, conf.level = conf.level, verbose = verbose,
                            search_pad_percent = 0.1,var.equal = var.equal)
  } else if (distribution =='z_dist') {
    ldm = ldm_normal_zdist(x,y, conf.level = conf.level, verbose = verbose,
                           search_pad_percent = 0.1,var.equal = var.equal)
  } else {
    ldm=NULL
    warning("Unrecognized distribution argument, only 't_dist', 'z_dist' allowed")
  }
  
  return(ldm)
   
}




ldm_normal_tdist <- function(x,y = NULL, conf.level = 0.95, verbose = FALSE, 
                             var.equal = FALSE, preserve_sign = FALSE) {
  #' @description  Calculate least difference in means statistic using a 
  #' t-distribution
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @return Returns Least Difference in Means
  #' @usage ldm_normal_tdist(x)
  #' @usage ldm_normal_tdist(x, y)
  #' @usage ldm_normal_tdist(x, y, conf.level = 0.95)
  #' @examples 
  #' x <- rnorm(n=8,mean=0,sd=1); 
  #' y <- rnorm(n=8,mean=0,sd=1); 
  #' ldm_normal_tdist(x,y)
  
  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distribution of the sampel for 1 sample)
  # n_x <- length(x); n_y <- length(y)
  # sd_x <- sd(x); sd_y <- sd(y)
  
  # Pooled mean, degrees of freedom, and standard deviation
  # if (is.null(y)) {
  #   # 1-sample distribution
  #   xbar_dm <- mean(x)
  #   sd_dm = sd_x / sqrt(n_x)
  # } else { 
  #   # 2-sample distribution
  #   xbar_dm <- mean(y) - mean(x)
  #   sd_pool = sqrt( ((nx-1)*sd_x^2 + (ny-1)*sd_y^2)/ (n_x+n_y-2) )
  #   sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
  # }
  
  # Calculate Confidence Limits
  t <- t.test(x=x,y=x, conf.level = conf.level, var.equal=var.equal)
  
  # Calculate ldm
  ldm = ldm_from_conf_int(t$conf.int[1], t$conf.int[2],preserve_sign = preserve_sign)
  
  # Calculate ldm
  return(ldm)  

}

ldm_normal_zdist <- function(x, y = NULL, conf.level = 0.95, verbose = FALSE, 
                             var.equal = FALSE, preserve_sign = TRUE) {
  #' @description Calculate least difference in meants using z-distribution
  #'
  #' @param x measurements from first group
  #' @param y measurements from second group (optional)
  #' @param conf.level confidence level for statistics (default 0.95)
  #' @param verbose function prints out extra input during execution
  #' @param var.equal boolean assume variance qual between groups (TODO: not 
  #' implemented for TRUE)
  #' @param preserve_sign flag to preserve sign of the least mean difference, 
  #' default TRUE
  #' @usage ldm_normal_zdist(x)
  #' @usage ldm_normal_zdist(x, y)
  #' @usage ldm_normal_zdist(x, y, conf.level = 0.95)
  #' @examples
  #' x <- rnorm(n=50,mean=0,sd=1); y <- rnorm(n=50,mean=0,sd=1); 
  #' ldm_normal_zdist(x,y)
  
  # Calculate basic stats of input samples defined by distribution d, the 
  
  # Calculate basic stats of input samples defined by distribution d, the 
  # difference distribution (or the distribution of the sample for 1 sample)
  n_x <- length(x); n_y <- length(y)
  sd_x <- sd(x); sd_y <- sd(y)
  alpha = 1- conf.level
  # Pooled mean, degrees of freedom, and standard deviation
  if (is.null(y)) {
    # 1-sample distribution
    xbar_dm <- mean(x)
    sd_dm = sd_x / sqrt(n_x)
  } else { 
    # 2-sample distribution
    xbar_dm <- mean(y) - mean(x)
    sd_pool = sqrt( ((n_x-1)*sd_x^2 + (n_y-1)*sd_y^2)/ (n_x+n_y-2) )
    sd_dm = sqrt( sd_x^2 / n_x  + sd_y^2 / n_y)
  }
  
  # Calculate z confidence intervals
  ci = qnorm(c(alpha/2, 1-alpha/2), mean = xbar_dm, sd= sd_dm)
  
  # Calculate ldm
  ldm = ldm_from_conf_int(ci[1], ci[2],preserve_sign = preserve_sign)
  
  return(ldm)  
}

ldm_from_conf_int <- function(ci_lower, ci_upper, preserve_sign = FALSE) {
  
  if (!preserve_sign) {
   ldm = (sign(ci_lower)==sign(ci_upper))*min(abs(c(ci_lower,ci_upper)))
  } else{
   ldm = sign(ci_lower) * (sign(ci_lower)==sign(ci_upper))*min(abs(c(ci_lower,ci_upper)))
  }
  
}