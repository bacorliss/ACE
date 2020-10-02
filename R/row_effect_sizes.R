

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


