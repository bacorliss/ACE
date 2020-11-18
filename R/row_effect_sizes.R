

# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)

p_load(BayesFactor)
p_load(TOSTER)


## User defined functions
# Calculate variance by row: sum of square of deviation from mean over n-1
rowVars <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)}
# Row standard deviation
rowSds  <- function(x, ...) sqrt(rowVars(x))
# Pooled standard deviation of two matrices, 1 sample per row
rowSdPooled <- function(m1,m2) sqrt( (rowSds(m1)^2 + rowSds(m2)^2 )/2 )
# Weighted pooled standard deviation of two matrices, 1 sample per row
rowSdPooledWeighted <- 
  function(m1,m2) sqrt( ( (dim(m1)[2] - 1) * rowVars(m1)   + 
                            (dim(m2)[2] - 1) * rowVars(m2) ) /
                          (dim(m1)[2] + dim(m2)[2] - 2    ) )


# Effect Size Statistics Functions
#
#------------------------------------------------------------------------------

# Delta Family Statistics
#
# Cohens D
#   d = (M2 - M1)/s_pool
row_cohend <- function (m1,m2) ( rowMeans(m1) - rowMeans(m2)) / rowSdPooled(m1,m2)

# Hedges G
#   d = (M2 - M1)/s_pool
row_hedgeg <- function (m1,m2) ( rowMeans(m1) - rowMeans(m2)) / rowSdPooledWeighted(m1,m2)

# Glass's Delta
#   d = (M2 - M1)/s_1
row_glassdelta <- function (m1,m2) ( rowMeans(m2) - rowMeans(m1)) / rowSds(m1) 


row_tscore_2s <- function(m1, m2) {
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  s1 <- rowSds(m1)
  s2 <- rowSds(m2)
  sm_pooled <-sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) *     (1/n1 + 1/n2) ) 
  tstat<-  ( rowMeans(m1) - rowMeans(m2)) / sm_pooled
}

row_zscore_2s <- function(m1, m2) {
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  s1 <- rowSds(m1)
  s2 <- rowSds(m2)
  zstat<-  ( rowMeans(m1) - rowMeans(m2)) /
    sqrt( s1^2/n1 + s2^2/n2)
}

row_ttest_2s <- function(m1, m2) {
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  s1 <- rowSds(m1)
  s2 <- rowSds(m2)
  sm_pooled <-sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) *     (1/n1 + 1/n2) ) 
  tstat<-  ( rowMeans(m1) - rowMeans(m2)) / sm_pooled
  p <- 2*pnorm(-abs(tstat))
}

row_ztest_2s <- function(m1, m2) {
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  s1 <- rowSds(m1)
  s2 <- rowSds(m2)
  zstat<-  ( rowMeans(m1) - rowMeans(m2)) /
    sqrt( s1^2/n1 + s2^2/n2)
  p  <- 2*pnorm(-abs(zstat))
}



row_mmd_2s_zdist <- function(m1, m2, ...) {
  mmd <- sapply(1:dim(m1)[1], function(i)  mmd_normal_zdist(m1[i,], m2[i,], ...))
}


row_bayesf_2s <- function(m1, m2, parallelize = FALSE, paired = FALSE) {
  # Ignore diagnostic messages during function call.
  wrap_fun <- function(x1,x2)   {
    suppressMessages(extractBF(ttestBF(x = x1,  y = x2, onlybf = TRUE))$bf)
  }
  # if (parallelize) {
  #   bf <- future_sapply(1:dim(m1)[1], function(i) wrap_fun(m1[i,], m2[i,]))
  # } else {
    bfx <- sapply(1:dim(m1)[1], function(i) wrap_fun(m1[i,], m2[i,]))
  # }
  # browser();
  return(bfx)
}


row_tost_2s_slow <- function (m1,m2) {
  
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  
  tost_fun <- function(x1,x2) 
    as.data.frame(dataTOSTtwo(
      data.frame(
        grp=as.factor(c(rep(1,n1),rep(2,n2))), 
        value=c(x1, x2)),
      deps="value", group = "grp", var_equal = FALSE, low_eqbound = -1e6,
      high_eqbound = 1e6, eqbound_type = "d", alpha = 0.05,
      desc = FALSE, plots = FALSE)$tost)$'p[0]'
  
  
  tost_p <- sapply(1:dim(m1)[1], function(i)   tost_fun(m1[i,], m2[i,]))
 
}


row_tost_2s <- function (m1,m2,low_eqbound = -1e-3,high_eqbound = 1e-3) {
  
  
  fast_2max <-function(v1,v2) {
    m = v1
    m[v2>v1] <- v2[v2>v1]
    return(m)
  }
  
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  s1 <- rowSds(m1)
  s2 <- rowSds(m2)
  s_md <- sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) * (1/n1 + 1/n2) )

  
  t_lower  <-  ( rowMeans(m1) - rowMeans(m2) -low_eqbound)  / s_md 
  p_low <- 1 - pt(t_lower, df = n1+n2-1)
  t_upper <-  ( rowMeans(m1) - rowMeans(m2) - high_eqbound)  / s_md 
  p_high <- pt(t_upper, df = n1+n2-1)
  # Quick max value of two vectors
  p_tost_fast <- fast_2max(p_low,p_high)
  p_tost_fast
  
  # # Test that mu[md] is greater than low
  # p_lower <- sapply(1:dim(m1)[1], function(i)
  #   t.test(m1[i,], m2[i,], alternative = "greater", mu = low_eqbound,
  #          paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value)
  # # Test that mu[md] lower than high
  # p_upper <- sapply(1:dim(m1)[1], function(i)
  #   t.test(m1[i,], m2[i,], alternative = "less", mu = high_eqbound,
  #          paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value)
  # p_tost_slow <- fast_2max( p_lower, p_upper)
  # p_tost_slow
  
  return(p_tost_fast)
}


## Test data
# 
# m1 = matrix(rnorm(1000, mean=0, sd=1), ncol = 50, nrow=20)
# m2 = matrix(rnorm(1000, mean=1, sd=1), ncol = 50, nrow=20)
# 
# # Z score and test
# z = rowzScore(m1, m2)
# p_row_ztest = 2*pnorm(-abs(z))
# p_actual_ztest = sapply(1:nrow(m1),function(x) t.test(m1[x,],m2[x,])$p.value)
# 
# # t score and test
# t = rowtScore(m1, m2)
# p_a = 2*pt(-abs(t), df = 98)
# 
# t_formal = sapply(1:nrow(m1),function(x) t.test(m1[x,],m2[x,])$stat)
# p_formal = sapply(1:nrow(m1),function(x) t.test(m1[x,],m2[x,])$p.value)




quantify_row_effect_sizes <- function(x_a, x_b, rand.seed = 0, 
                                  parallelize_bf = FALSE) {
  set.seed(rand.seed)
  
  # Use Exp 1 and 2 coefficients to generate data from normalized base data
  x_a = matrix(rnorm(df$n_samples * df$n_1a, mean = df$mu_1a, 
                     sd = df$sigma_1a), nrow = df$n_samples, 
               ncol = df$n_1a)
  x_b = matrix(rnorm(df$n_samples * df$n_1b, mean = df$mu_1b, 
                     sd = df$sigma_1b), nrow = df$n_samples, 
               ncol = df$n_1b)
  
  # Means
  xdbar_d = rowMeans(x_b) - rowMeans(x_a)
  df$mean_xdbar_d = mean(xdbar_d)
  df$sd_xdbar = sd(xdbar_d)
  
  # Stds
  s_md = sqrt(rowSds(x_a)^2/df$n_a + rowSds(x_b)^2/df$n_b)
  df$mean_sdmd = mean(s_md)
  df$sd_sdmd   = sd(s_md)
  
  # Rel Means: mean divided by control mean
  rxdbar = xdbar_d/rowMeans(x_a)
  df$mean_rxdbar = mean(rxdbar)
  df$sd_rxdbar   = sd(rxdbar)
  
  # Rel STDs: sd divided by control mean
  # Relaltive STD is halfway between xbar_A and xbar_B
  exp1_rsd_md = s_1md / (rowMeans(x_a) + 0.5 * xdbar_1d)
  df$exp1_mean_rsdmd = mean(exp1_rsd_md)
  df$exp1_sd_rsdmd   = sd(exp1_rsd_md)
  
  # Bayes Factor
  diff_bf = row_bayesf_2s(x_a, x_b, parallelize = parallelize_bf, paired = FALSE)
  df$exp_mean_bf = mean(diff_bf)
  
  # NHST P-value 
  # The more equal experiment will have a larger p-value
  diff_z_score <- row_zscore_2s(x_b, x_a)
  diff_pvalue = 2*pnorm(-abs(diff_z_score))
  df$exp_mean_pvalue = mean(diff_pvalue)
  df$exp_sd_pvalue = sd(diff_pvalue)
  
  # TOST p value (Two tailed equivalence test)
  diff_tostp <- row_tost_2s(x_b, x_a,low_eqbound = -.1,high_eqbound = .1)
  df$exp_mean_tostp = mean(diff_tostp)
  df$exp_sd_tostp = sd(diff_tostp)
  
  
  # Cohens D
  diff_cohend = row_cohend(x_a, x_b)
  df$exp1_mean_cohend = mean(diff_cohend)
  df$exp1_sd_cohend = sd(diff_cohend)
  
  # Most Mean Diff
  diff_mmd = row_mmd_2s_zdist(x_a, x_b)
  df$exp_mean_mmd = mean(diff_mmd)
  df$exp_sd_mmd = sd(diff_mmd)
  
  # Relative Most Mean Diff
  diff_rmmd = diff_mmd / rowMeans(x_a)
  df$exp_mean_rmmd = mean(diff_rmmd)
  df$exp_sd_rmmd = sd(diff_rmmd)
  
  # Random group
  x_rand = rowMeans(matrix(rnorm(df$n_samples * df$df_1d, mean = 0, sd = 1), 
                           nrow = df$n_samples, ncol = df$df_1d))
  df$exp_mean_x_rand = mean(x_rand)
  df$exp_sd_x_rand = sd(x_rand)
  
  
  return(df)
}

