


#' Library of functions to compute various 2-sample effect size statistics in a 
#' row-by-row fashion. Input are two mxn matrices with m samples as rows and n 
#' measurements as columns as columns. All functions return a vector of length m, 
#' with an effect size quantified for each sample.

# Test input for functions
# m1 = matrix(rnorm(1000, mean=10, sd=1), ncol = 50, nrow=20)
# m2 = matrix(rnorm(1000, mean=15, sd=1), ncol = 50, nrow=20)

# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(BayesFactor)
p_load(TOSTER)
# Second generation p-values
p_load(sgpv)

source("R/ldm.R")
source("R/mdm.R")
source("R/rationormal_toolbox.R")

## User defined functions


# Find min or max value between two columns: optimized for speed since its a 
# smaller edge case than the typical rowMin/rowMax application
rowmin_2col <- function(v1,v2) { 
  if(length(v2)==1) { v1[v2<v1] <- v2
  } else { v1[v2<v1] <- v2[v2<v1]}
  return(v1)}
rowmax_2col <- function(v1,v2) { v1[v2>v1] <- v2[v2>v1]; return(v1)}

# Calculate variance by row: sum of square of deviation from mean over n-1
rowVars <- function(x, ...) {rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1)}
# Row standard deviation
rowSds  <- function(x, ...) sqrt(rowVars(x))
# Pooled standard deviation of two matrices, 1 sample per row
rowSdCombined <- function(m1,m2) sqrt( (rowSds(m1)^2 + rowSds(m2)^2 )/2 )
# Weighted pooled standard deviation of two matrices, 1 sample per row
rowSdPooled <- 
  function(m1,m2) 
    sqrt( ( (dim(m1)[2] - 1) * rowVars(m1) +  
            (dim(m2)[2] - 1) * rowVars(m2) ) /
                          (dim(m1)[2] + dim(m2)[2] - 2    ) )


# Delta Family Statistics
#
# Cohens D
#   d = (M2 - M1)/s_pool
row_cohend <- function (m1,m2) ( rowMeans(m1) - rowMeans(m2)) / rowSdCombined(m1,m2)

# Hedges G
#   d = (M2 - M1)/s_pool
row_hedgeg <- function (m1,m2) ( rowMeans(m1) - rowMeans(m2)) / rowSdPooled(m1,m2)

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


row_confint_z  <- function(m1,m2, conf.level = 0.95) {
  x_bar_1 = rowMeans(m1); x_bar_2 = rowMeans(m2)
  s_pool = rowSdPooled(m1, m2)
  x_dm <- x_bar_2 - x_bar_1
  
  df_ci <- data.frame(
    conf.lo = x_dm - qnorm(1 - (1-conf.level)/2) * s_pool * sqrt(1/dim(m1)[2] + 1/dim(m2)[2]),
    conf.hi = x_dm + qnorm(1 - (1-conf.level)/2) * s_pool * sqrt(1/dim(m1)[2] + 1/dim(m2)[2]) )
  return(df_ci)
}

row_confint_t  <- function(m1,m2, conf.level = 0.95) {
  ci <- sapply(1:dim(m1)[1], function(i)  t.test(x=m2[i,], y=m1[i,])$conf.int)
  df_ci <- data.frame(conf.lo = ci[1,], conf.hi = ci[2,])
  return(df_ci)
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


row_ci_mean_2s_zdist <-function(m1,m2, conf.level=0.95) {
  #' @description Calculates row by row z distribution confidence interval of 
  #' the mean
  #' @param m1 matrix of measurements, with rows as samples and columns as 
  #' measurements
  #' @param m2 matrix of measurements, with rows as samples and columns as 
  #' measurements
  #' @param conf.level confidence level for confidence interval
  #' @return 

  xbar_dm <- rowMeans(m2) - rowMeans(m1)
  z_stat <- qnorm(1-(1-conf.level)/2)
  x_offset <- z_stat * rowSdPooled(m1,m2) * sqrt( 1/dim(m1)[2] + 1/dim(m2)[2] )
  
  ci = tibble(ci_lower = xbar_dm - x_offset, ci_upper = xbar_dm + x_offset)
  return(ci)
}

row_mdm_2s_zdist <- function(m1, m2, ...) {
  mdm <- sapply(1:dim(m1)[1], function(i)  mdm_normal_zdist(m1[i,], m2[i,], ...))
}




row_rmdm_2s_zdist <- function(m1, m2, mdms = NULL, conf.level = 0.95, method = "qnormrat") {
  #' @description calculates relative mdm row by row given a matrix of control 
  #' samples and experiment samples (limitation: samples must have same sample 
  #' size within groups). M1 and M2 must have same number of rows (samples), but 
  #' can have different numbers of columns (measurements)
  #'
  #' @param m1 control group
  #' @param m2 experiment group
  #' 
  #' @return vector of rmdm values, one for each row of m1 and m2
  # browser()

  # Calculate mdm from data if not supplied
  # if (is.null(mdms)) {
  #    mdms <- sapply(1:dim(m1)[1], function(i)  mdm_normal_zdist(m1[i,], m2[i,], conf.level = conf.level.mdm))
  # }
  # Calculate standard deviation of difference in means
  # s_dm <- sqrt(rowvars(m1)/dim(m1)[2] + rowvars(m2)/dim(m2)[2])
  # # Calculate means and standard error of control group
  # # control
  # xbar_1 <- rowmeans(m1)
  # se_1 <- rowsds(m1)/dim(m1)[2]

  if (method =="qnormrat") {
    # Calculate mean and std of difference in means distribution
    
    
    source("R/rationormal_toolbox.R")
    rmdms <- sapply(1:dim(m1)[1], function(i)  
      rmdm_normal_zdist(x = m1[i,], y = m2[i,], mdm = NULL, conf.level, method = "qnormrat"))

   
  } else if (method =="mdm_normalized") {
    
    mdm <- sapply(1:dim(m1)[1], function(i)  
      mdm_normal_zdist(x=m1[i,], y=m2[i,], conf.level = sqrt(0.05)))
    rmdms <- sapply(1:dim(m1)[1], function(i)  
      rmdm_normal_zdist(x=m1[i,], y=m2[i,], mdm = mdms[i], conf.level), method = method)
    
  } else {stop('row_rmdm_2s_zdist(): unsupported method')}
  
  

  
  return(rmdms)
}


row_rmdm_normal_montecarlo <- function(mc, me, conf.level = 0.95, n_trials=1e6) {
  
  
  n_samples = dim(mc)[1]
    
  means_c = rowMeans(mc)
  n_c = dim(mc)[2]
  sds_c = rowSds(mc)
  ses_c = sds_c/sqrt(n_c)
  
  means_e = rowMeans(me)
  n_e = dim(me)[2]
  sds_e = rowSds(me)
  ses_e = sds_e/sqrt(n_e)
  
  means_dm = means_e - means_c
  ses_dm = sqrt(sds_c^2/n_c + sds_e^2/n_e)
  
  
  # Generate seed random numbers to save time between rows
  seed_1 = rnorm(n_trials, mean = 0,  sd = 1)
  seed_2 = rnorm(n_trials, mean = 0,  sd = 1)
  
  
  # time1 <- Sys.time()
  # means_c_mat = seed_1 %*% t(sds_c) + rep(means_c, each = n_samples)
  # means_e_mat = seed_2 %*% t(sds_e) + rep(means_e, each = n_samples)
  # 
  # rxbar_mat = abs(means_e_mat - means_c_mat)/means_c_mat
  # # Compute monte carlo in matrix format to avoid for loop
  # 
  # # Calculate quantile for each column
  # rmdms = apply(rxbar_mat , 2 , quantile, probs = conf.level , na.rm = TRUE )
  # time2 <- Sys.time()
  # time2 - time1
  
  
  # time1 <- Sys.time()
  rmdms <- sapply(1:dim(mc)[1], function(i)  
    rmdm_normal_montecarlo(means_dm[i], ses_dm[i], means_c[i], ses_c[i], 
                           abs(seed_1 * ses_e[i] + means_e[i] - seed_2  * ses_c[i] + means_c[i]), 
                           seed_2  * ses_c[i] + means_c[i],
                           conf.level = 0.95, n_trials=n_trials))
  # time2 <- Sys.time()
  # time2 - time1
  # 
  return(rmdms)
}



row_ratio_normal_coe <- function(m1, m2, conf.level = 0.95, method = "qnormrat") {
  #' @description calculates relative mdm row by row given a matrix of control 
  #' samples and experiment samples (limitation: samples must have same sample 
  #' size within groups). M1 and M2 must have same number of rows (samples), but 
  #' can have different numbers of columns (measurements)
  #'
  #' @param m1 control group
  #' @param m2 experiment group
  #' 
  #' @return vector of rmdm values, one for each row of m1 and m2
  # browser()
  
  # Calcualt mdm from data if not supplied
  # if (is.null(mdms)) {
  #    mdms <- sapply(1:dim(m1)[1], function(i)  mdm_normal_zdist(m1[i,], m2[i,], conf.level = conf.level.mdm))
  # }
  # Calculate standard deviation of difference in means
  # s_dm <- sqrt(rowvars(m1)/dim(m1)[2] + rowvars(m2)/dim(m2)[2])
  # # Calculate means and standard error of control group
  # # control
  # xbar_1 <- rowmeans(m1)
  # se_1 <- rowsds(m1)/dim(m1)[2]
  
  # Numerator: experiment sample
  means_x <- rowMeans(m2)
  n_x <- dim(m2)[2]
  sds_x <- rowSds(m2)
  ses_x <- sds_x/sqrt(n_x)
  # sds_x <- rowSds(m1)/sqrt(n_x)
  
  # Denominator: control sample
  means_y <- rowMeans(m1)
  n_y <- dim(m1)[2]
  sds_y <- rowSds(m1)
  ses_y = sds_y/sqrt(n_y)
  # sds_y <- rowSds(m2)/ sqrt(n_y)
  
  xbar_dm = means_x - means_y
  sd_dm = sqrt(sds_x^2/n_x + sds_y^2/n_y)
  df_dm = n_x + n_y - 2
   
  # browser();
  
  if (method == 'ttestratio_default') {
    ratio_cls <- 
      abs(t(sapply(1:dim(m1)[1], function(i) 
        ttestratio_default(x=m2[i,], y=m1[i,],alternative = "two.sided", rho = 1, 
                           var.equal = TRUE, conf.level = 1 - 2*(1-conf.level))$conf.int)))
    ratio_ucl <- rowmax_2col(ratio_cls[,1],ratio_cls[,2])
    
  } else if (method == 'ttestratio') {
    ratio_cls <- t(sapply(1:dim(m1)[1], function(i)
      ttestratio(mx = means_x[i], sdx = sds_x[i], dfx = n_x-1,
                 my = means_y[i], sdy = sds_y[i], dfy = n_y-1,
                 alternative = "two.sided", rho = 1, var.equal = TRUE,
                 conf.level =  1 - 2*(1-conf.level))$conf.int))
    ratio_ucl <- rowmax_2col(ratio_cls[,1],ratio_cls[,2])
    
  } else if (method == 'qnormrat') {
    
    browser();
    
    # Works with positive values mu_b/mu_a
    ratio_cl_lo <- sapply(1:dim(m1)[1], function(i)
      qnormrat(p = 1-conf.level, means_x[i], ses_x[i], means_y[i], ses_y[i],
               VERBOSE=FALSE))
    
    ratio_cl_hi <- sapply(1:dim(m1)[1], function(i)
      qnormrat(p = conf.level, means_x[i], ses_x[i], means_y[i], ses_y[i],
               VERBOSE=FALSE))
    
    ratio_ucl <- rowmax_2col(abs(ratio_cl_lo), abs(ratio_cl_hi))

  } else {stop('row_ratio_normal_coe(): unsupported method')}
  
  
  
  return(ratio_ucl)
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




row_ldm_2s_zdist <- function(m1, m2,...) {
  ldm <- sapply(1:dim(m1)[1], function(i)  ldm_normal_zdist(x=m1[i,], y=m2[i,], ...))
  return(ldm)
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


row_sgpv <- function(m1, m2, null.lo, null.hi){
  df_ci <- row_confint_t(m1, m2)
  p_sg <- sgpvalue( df_ci$conf.lo, df_ci$conf.hi, null.lo = null.lo, null.hi = null.hi)$p.delta
  return(p_sg)
}



row_tost_2s <- function (m1,m2,low_eqbound = -1e-3,high_eqbound = 1e-3, conf.level = 0.95) {
  
  n1 <- dim(m1)[2]
  n2 <- dim(m2)[2] 
  s1 <- rowSds(m1)
  s2 <- rowSds(m2)
  s_dm <- sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) * (1/n1 + 1/n2) )

  
  t_lower  <-  ( rowMeans(m1) - rowMeans(m2) -low_eqbound)  / s_dm 
  p_low <- rowmin_2col((1 - pt(t_lower, df = n1+n2-1)) * 0.05/(1-conf.level),1)
  
  t_upper <-  ( rowMeans(m1) - rowMeans(m2) - high_eqbound)  / s_dm 
  p_high <- rowmin_2col(pt(t_upper, df = n1+n2-1) * 0.05/(1-conf.level),1)
  # Quick max value of two vectors
  p_tost_fast <- rowmax_2col(p_low,p_high)
  
  
  # # Test that mu[md] is greater than low
  # p_lower <- sapply(1:dim(m1)[1], function(i)
  #   t.test(m1[i,], m2[i,], alternative = "greater", mu = low_eqbound,
  #          paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value)
  # # Test that mu[md] lower than high
  # p_upper <- sapply(1:dim(m1)[1], function(i)
  #   t.test(m1[i,], m2[i,], alternative = "less", mu = high_eqbound,
  #          paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value)
  # p_tost_slow <- rowmax_2col( p_lower, p_upper)
  
  return(p_tost_fast)
}


## Test data
# # 
m1 = matrix(rnorm(1000, mean=10, sd=1), ncol = 50, nrow=20)
m2 = matrix(rnorm(1000, mean=15, sd=1), ncol = 50, nrow=20)
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




quantify_row_stats <- function(x_a, x_b, parallelize_bf = FALSE, stat_exclude_list=NULL, conf.level = 0.95) {
  #' @description Given two matrices of measurements, with rows representing 
  #' samples and columns observations, calculates a collection of effect size
  #' statistics and reports the mean and standard deviation across samples (only 
  #' supports same sample size within x_a and x_b). Sign is preserved with the 
  #' statistics, when comparing them for agreement and disagreement, the 
  #' magnitude is taken (parent function).
  #' 
  #' @param x_a matrix of observations from control group where rows are separate
  #'  samples
  #' @param x_b matrix of observations from experiment group where rows are 
  #' separate samples 
  #' @param parallelize_bf flag for parallel processing of Bayes Factor 
  #' (since current implementation in R is slow)
  #' @param exclude_list list of variables to exclude
  #' 
  #' @return dataframe of mean and standard deviation for each effect size statistic
  
  n_a = dim(x_a)[2]
  n_b = dim(x_b)[2]
  n_samples = dim(x_a)[1]
  df_d = n_a + n_b - 2
  
  # Initialize data frames so that
  stat_list <- c("xbar_dm", "rxbar_dm", "sd_dm", "rsd_dm", "bf", "pvalue", "tostp",
                 "cohend", "mdm", "rmdm", "ldm", "rldm", "rand")
  
  df = data.frame(matrix(ncol = length(stat_list), nrow=n_samples))
  colnames(df) <- stat_list
  df_hat <- df[1,]
  df_hdt <- df[1,]
  df_name <- df[1,]
  df_pretty <- df[1,]
    
  # 1) Mean of the difference of means
  df$xbar_dm = rowMeans(x_b) - rowMeans(x_a)
  df_hat$xbar_dm <- "<"
  df_hdt$xbar_dm <- ">"
  df_pretty$xbar_dm <- "bar(x)[DM]"

  # 2) Rel Means: mean of the difference in means divided by control group mean
  df$rxbar_dm = df$xbar_dm/rowMeans(x_a)
  df_hat$rxbar_dm <- "<"
  df_hdt$rxbar_dm <- ">"
  df_pretty$rxbar_dm <- "r*bar(x)[DM]"
  
  # 3) Std of the difference in means
  df$sd_dm = sqrt(rowSds(x_a)^2/n_a + rowSds(x_b)^2/n_b)
  df_hat$sd_dm <- "<"
  df_hdt$sd_dm <- "<"
  df_pretty$sd_dm <- "s[DM]"
  
  # 4) Relative STD: std of difference in means divided by midpoint between group means
  df$rsd_dm = df$sd_dm / (rowMeans(x_a) + 0.5 * df$xbar_dm)
  df_hat$rsd_dm <- "<"
  df_hdt$rsd_dm <- "<"
  df_pretty$rsd_dm <- "r*s[DM]"
  
  # 5) Bayes Factor
  df$bf = row_bayesf_2s(x_a, x_b, parallelize = parallelize_bf, paired = FALSE)
  df_hat$bf <- "<"
  df_hdt$bf <- ">"
  df_pretty$bf <- "BF"
  
  # 6) NHST P-value: The more equal experiment will have a larger p-value
  diff_z_score <- row_zscore_2s(x_b, x_a)
  df$pvalue = rowmin_2col(v1 = 2*pnorm(-abs(diff_z_score))*0.05/(1-conf.level),1)
  df_hat$pvalue <- ">"
  df_hdt$pvalue <- "<"
  df_pretty$pvalue <- "p[N]"
  
  # 7) TOST p value (Two tailed equivalence test)
  df$tostp <- row_tost_2s(x_b, x_a,low_eqbound = -.1,high_eqbound = .1, conf.level = conf.level)
  df_hat$tostp <- "<"
  df_hdt$tostp <- ">"
  df_pretty$tostp <- "p[E]"
  
  # 8) Cohens D
  df$cohend = row_cohend(x_a, x_b)
  df_hat$cohend <- "<"
  df_hdt$cohend <- ">"
  df_pretty$cohend <- "Cd"
  
  # 9) most difference in means
  df$mdm = row_mdm_2s_zdist(x_a, x_b, conf.level = conf.level)
  df_hat$mdm <- "<"
  df_hdt$mdm <- ">"
  df_pretty$mdm <- "delta[M]"
  
  # 10) Relative most difference in means
  df$rmdm = df$mdm / rowMeans(x_a)
  df_hat$rmdm <- "<"
  df_hdt$rmdm <- ">"
  df_pretty$rmdm <- "r*delta[M]"
  
  # 11) Least Difference in Means
  df$ldm = row_ldm_2s_zdist(x_a, x_b, conf.level=conf.level)
  df_hat$ldm <- "<"
  df_hdt$ldm <- ">"
  df_pretty$ldm <- "delta[L]"
  
  # 12) Relative Least Difference in Means
  df$rldm = df$ldm / rowMeans(x_a)
  df_hat$rldm <- "<"
  df_hdt$rldm <- ">"
  df_pretty$rldm <- "r*delta[L]"
  
  # 13) Random group
  df$rand = rowMeans(matrix(rnorm(n_samples * 50, mean = 0, sd = 1), 
                           nrow = n_samples, ncol = 50))
  df_hat$rand <- "<"
  df_hdt$rand <- ">"
  df_pretty$rand <- "Rnd"
  
  
  # Drop any statistics requested by user
  df <- df[ , !(names(df) %in% stat_exclude_list)]
  df_hat <- df_hat[ , !(names(df_hat) %in% stat_exclude_list)]
  df_hdt <- df_hdt[ , !(names(df_hdt) %in% stat_exclude_list)]
  df_pretty <- df_pretty[ , !(names(df_pretty) %in% stat_exclude_list)]
  
  # Store attributes within df
  attr(df,"user_attributes") <- c("user_attributes","hat", "varnames", "varnames_pretty")
  attr(df,"hat") <- df_hat
  attr(df,"hdt") <- df_hdt
  attr(df,"varnames") <- colnames(df)
  attr(df,"varnames_pretty") <- df_pretty
  
  
  # dfs = list("df" = df, "df_hat" = df_hat)
  return(df)
}

