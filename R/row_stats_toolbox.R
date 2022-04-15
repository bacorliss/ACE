


#' Library of functions to compute various 2-sample effect size statistics in a 
#' row-by-row fashion. Input are two mxn matrices with m samples as rows and n 
#' measurements as columns as columns. All functions return a vector of length m, 
#' with an effect size quantified for each sample.

# Test input for functions
# m_c = matrix(rnorm(1000, mean=10, sd=1), ncol = 50, nrow=20)
# m_e = matrix(rnorm(1000, mean=15, sd=1), ncol = 50, nrow=20)

# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(BayesFactor)
p_load(TOSTER)
# Second generation p-values
p_load(sgpv)
source("R/aces.R")

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
rowSdCombined <- function(m_c,m_e) sqrt( (rowSds(m_c)^2 + rowSds(m_e)^2 )/2 )
# Weighted pooled standard deviation of two matrices, 1 sample per row
rowSdPooled <- 
  function(m_c,m_e) 
    sqrt( ( (dim(m_c)[2] - 1) * rowVars(m_c) +  
            (dim(m_e)[2] - 1) * rowVars(m_e) ) /
                          (dim(m_c)[2] + dim(m_e)[2] - 2    ) )


# Delta Family Statistics
#
# Cohens D
#   d = (m_e - m_c)/s_pool
row_cohend <- function (m_c,m_e) ( rowMeans(m_c) - rowMeans(m_e)) / rowSdCombined(m_c,m_e)

# Hedges G
#   d = (m_e - m_c)/s_pool
row_hedgeg <- function (m_c,m_e) ( rowMeans(m_c) - rowMeans(m_e)) / rowSdPooled(m_c,m_e)

# Glass's Delta
#   d = (m_e - m_c)/s_1
row_glassdelta <- function (m_c,m_e) ( rowMeans(m_e) - rowMeans(m_c)) / rowSds(m_c) 


row_tscore_2s <- function(m_c, m_e) {
  n1 <- dim(m_c)[2]
  n2 <- dim(m_e)[2] 
  s1 <- rowSds(m_c)
  s2 <- rowSds(m_e)
  sm_pooled <-sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) *     (1/n1 + 1/n2) ) 
  tstat<-  ( rowMeans(m_c) - rowMeans(m_e)) / sm_pooled
}

row_zscore_2s <- function(m_c, m_e) {
  n1 <- dim(m_c)[2]
  n2 <- dim(m_e)[2] 
  s1 <- rowSds(m_c)
  s2 <- rowSds(m_e)
  zstat<-  ( rowMeans(m_c) - rowMeans(m_e)) /
    sqrt( s1^2/n1 + s2^2/n2)
}


row_confint_z  <- function(m_c,m_e, conf.level = 0.95) {
  x_bar_1 = rowMeans(m_c); x_bar_2 = rowMeans(m_e)
  s_pool = rowSdPooled(m_c, m_e)
  x_dm <- x_bar_2 - x_bar_1
  
  df_ci <- data.frame(
    conf.lo = x_dm - qnorm(1 - (1-conf.level)/2) * s_pool * sqrt(1/dim(m_c)[2] + 1/dim(m_e)[2]),
    conf.hi = x_dm + qnorm(1 - (1-conf.level)/2) * s_pool * sqrt(1/dim(m_c)[2] + 1/dim(m_e)[2]) )
  return(df_ci)
}

row_confint_t  <- function(m_c,m_e, conf.level = 0.95) {
  ci <- sapply(1:dim(m_c)[1], function(i)  t.test(x=m_e[i,], y=m_c[i,], conf.level = conf.level) $conf.int)
  df_ci <- data.frame(conf.lo = ci[1,], conf.hi = ci[2,])
  return(df_ci)
}



row_ttest_2s <- function(m_c, m_e, conf.level) {
  p <- sapply(1:dim(m_c)[1], function(i)  
    t.test(x=m_e[i,],y=m_c[i,], conf.level =conf.level, var.equal = FALSE)$p.value)
}

row_ztest_2s <- function(m_c, m_e) {
  n1 <- dim(m_c)[2]
  n2 <- dim(m_e)[2] 
  s1 <- rowSds(m_c)
  s2 <- rowSds(m_e)
  zstat<-  ( rowMeans(m_c) - rowMeans(m_e)) /
    sqrt( s1^2/n1 + s2^2/n2)
  p  <- 2*pnorm(-abs(zstat))
}


row_ci_mean_2s_zdist <-function(m_c,m_e, conf.level=0.95) {
  #' @description Calculates row by row z distribution confidence interval of 
  #' the mean
  #' @param m_c matrix of measurements, with rows as samples and columns as 
  #' measurements
  #' @param m_e matrix of measurements, with rows as samples and columns as 
  #' measurements
  #' @param conf.level confidence level for confidence interval
  #' @return 

  xbar_dm <- rowMeans(m_e) - rowMeans(m_c)
  z_stat <- qnorm(1-(1-conf.level)/2)
  x_offset <- z_stat * rowSdPooled(m_c,m_e) * sqrt( 1/dim(m_c)[2] + 1/dim(m_e)[2] )
  
  ci = tibble(ci_lower = xbar_dm - x_offset, ci_upper = xbar_dm + x_offset)
  return(ci)
}

row_mdm <- function(m_c, m_e, conf.level, rand.seed = NA) {
  mdm <- sapply(1:dim(m_c)[1], function(i)  mdm_credint(m_e[i,], m_c[i,],
          conf.level = conf.level, relative = FALSE, rand.seed = rand.seed))

  return(mdm)
}


row_rmdm <- function(m_c, m_e, mdms = NULL, conf.level = 0.95, rand.seed = NA) {
  #' @description calculates relative mdm row by row given a matrix of control 
  #' samples and experiment samples (limitation: samples must have same sample 
  #' size within groups). m_c and m_e must have same number of rows (samples), but 
  #' can have different numbers of columns (measurements)
  #'
  #' @param m_c control group
  #' @param m_e experiment group
  #' 
  #' @return vector of rmdm values, one for each row of m_c and m_e
    rmdms <-  sapply(1:dim(m_c)[1], function(i)
      mdm_credint(x = m_e[i,], y = m_c[i,], conf.level = conf.level,
                  relative = TRUE, rand.seed = rand.seed))
  return(rmdms)
}

row_ldm <- function(m_c, m_e, conf.level, rand.seed = NA) {
  # browser()
  ldm <- sapply(1:dim(m_c)[1], function(i)  
    ldm_credint(m_e[i,], m_c[i,], conf.level = conf.level, relative = FALSE, 
                rand.seed = rand.seed))
  return(ldm)
}

row_rldm <- function(m_c, m_e, conf.level = 0.95, rand.seed = NA) {
  # browser()
  rldms <-  sapply(1:dim(m_c)[1], function(i)
    ldm_credint(x = m_e[i,], y = m_c[i,], conf.level = conf.level, 
                relative = TRUE, rand.seed = rand.seed))
  return(rldms)
}


row_macb_tdist_2sample  <- function(m_c, m_e, ...) {
  macb <- sapply(1:dim(m_c)[1], function(i)  macb_tdist_2sample (m_c[i,], m_e[i,], ...))
  return(macb)
}


quiet <- function(x) { sink(tempfile()); on.exit(sink()); invisible(force(x)) } 

row_bayesf_2s <- function(m_c, m_e, deltas = NULL, paired = FALSE) {
  #' @description Calculates bayes factor across rows for two groups.
  #' @param m_c matrix of observations from control group where rows are separate
  #'  samples
  #' @param m_e matrix of observations from experiment group where rows are 
  #' separate samples 
  #' @param nullInterval null region bound for Bayes factor in Standardized units
  #' @param paired boolean for paired data
  #' 
  #' @return vector of bayes factor 1xN samples
  #' 
  # Ignore diagnostic messages during function call.
  
  # If only one delta specified, make it delta for all rows
  if (length(deltas)==1) {deltas <- rep(deltas, dim(m_c)[2])}
  
  wrap_fun <-  function(x1,x2, nullInterval)   
    {extractBF(ttestBF(x = x1,  y = x2, nullInterval = nullInterval, paired  = paired), onlybf = TRUE)}[1]
  capture.output(type="message",
                 bfx <- sapply(1:dim(m_c)[1], function(i) 
                   wrap_fun(m_c[i,], m_e[i,], c(-deltas[i],deltas[i]))))
  
  return(bfx)
}



# row_tost_2s_slow <- function (m_c,m_e, low_eqbound, high_eqbound, conf.level) {
#   
#   n1 <- dim(m_c)[2]
#   n2 <- dim(m_e)[2] 
#   
#   tost_fun <- function(x1,x2) 
#     as.data.frame(dataTOSTtwo(
#       data.frame(
#         grp=as.factor(c(rep(1,n1),rep(2,n2))), 
#         value=c(x1, x2)),
#       deps="value", group = "grp", var_equal = FALSE, low_eqbound = -1e6,
#       high_eqbound = 1e6, eqbound_type = "d", alpha = 1-conf.level,
#       desc = FALSE, plots = FALSE)$tost)$'p[0]'
#   
#   tost_p <- sapply(1:dim(m_c)[1], function(i)   tost_fun(m_c[i,], m_e[i,]))
#  
# }


row_sgpv <- function(m_c, m_e,  null.los, null.his, conf.level){
  #' @description Calculates second generation p-values across rows for two groups.
  #' @param m_c matrix of observations from control group where rows are separate
  #'  samples
  #' @param m_e matrix of observations from experiment group where rows are 
  #' separate samples 
  #' @param null.lo null region bound for Bayes factor in Standardized units
  #' @param null.hi list of variables to exclude
  #' 
  #' @return vector of p-values, 1xN, where N is number of samples
  #' 

  df_ci <- row_confint_t(m_c, m_e, conf.level = conf.level)
  
  p_sg <- sgpvalue( df_ci$conf.lo, df_ci$conf.hi, null.lo = null.los,
                    null.hi = null.his)$p.delta
  return(p_sg)
}



row_tost_2s <- function (m_c, m_e, deltas = NULL, conf.level = 0.95) {

  # browser()
  n1 <- dim(m_c)[2]
  n2 <- dim(m_e)[2] 
  s1 <- rowSds(m_c)
  s2 <- rowSds(m_e)
  s_dm <- sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) * (1/n1 + 1/n2) )

  # If only one delta specified, make it delta for all rows
  if (length(deltas)==1) {deltas <- rep(deltas, dim(m_c)[1])}
  # Deltas define lower and upper bound for each sample
  lo_bounds <- -deltas
  hi_bounds <- deltas
  
  t_lower  <-  ( rowMeans(m_c) - rowMeans(m_e) - lo_bounds)  / s_dm 
  p_low <- rowmin_2col((1 - pt(t_lower, df = n1+n2-1)) * 0.05/(1-conf.level),1)
  
  t_upper <-  ( rowMeans(m_c) - rowMeans(m_e) - hi_bounds)  / s_dm 
  p_high <- rowmin_2col(pt(t_upper, df = n1+n2-1) * 0.05/(1-conf.level),1)
  # Quick max value of two vectors
  p_tost_fast <- rowmax_2col(p_low,p_high)
  
  
  # # Test that mu[md] is greater than low
  # p_lower <- sapply(1:dim(m_c)[1], function(i)
  #   t.test(m_c[i,], m_e[i,], alternative = "greater", mu = low_eqbound,
  #          paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value)
  # # Test that mu[md] lower than high
  # p_upper <- sapply(1:dim(m_c)[1], function(i)
  #   t.test(m_c[i,], m_e[i,], alternative = "less", mu = high_eqbound,
  #          paired = FALSE, var.equal = FALSE, conf.level = 0.95)$p.value)
  # p_tost_slow <- rowmax_2col( p_lower, p_upper)
  
  
  # pe = TOSTtwo.raw(mean(xa), mean(xb), sd(xa), sd(xb), length(xa), length(xb), 
  #             low_eqbound = -1e-3, high_eqbound = 1e-3, alpha = parse_fract(df[[alphadm_str]][n]),
  #             var.equal = FALSE, plot = FALSE, verbose = FALSE)
  
  
  return(p_tost_fast)
}


## Test data
# m_c = matrix(rnorm(1000, mean=10, sd=1), ncol = 50, nrow=20)
# m_e = matrix(rnorm(1000, mean=15, sd=1), ncol = 50, nrow=20)




quantify_row_stats <- function(x_a, x_b, parallelize_bf = FALSE, 
                               stat_exclude_list=NULL, conf.level = 0.95, delta,
                               is_delta_relative) {
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
  #' @param delta the half width of the null region used for many of the 
  #' statistics, including tostp, p2, and bf. Null region assumed to be zero
  #' centered, so a null region of [-1,1] has a delta of 1.Units assumed to be in
  #' units of effect size, so no normalized to pooled standard deviation.
  #' 
  #' @return dataframe of mean and standard deviation for each effect size statistic
  
  # save(list = ls(all.names = TRUE), file = "temp/quantify_row_stats.RData",envir = environment())
  # load(file = "temp/quantify_row_stats.RData")
  
  n_a = dim(x_a)[2]
  n_b = dim(x_b)[2]
  n_samples = dim(x_a)[1]
  df_d = n_a + n_b - 2

  # Deltas defines distance from zero in unscaled or relative units
  if (is_delta_relative) {deltas = delta*rowMeans(x_a)} else { deltas = rep(delta, n_samples)}
  
  # Initialize data frames so that only included stats are columns
  stat_list <- setdiff(c("xbar_dm", "rxbar_dm", "sd_dm", "rsd_dm", "bf", "pvalue",
                         "tostp", "p2", "cohend", "mdm", "rmdm", "ldm", "rldm",
                         "rand"), stat_exclude_list)
  df = data.frame(matrix(ncol = length(stat_list), nrow=n_samples))
  colnames(df) <- stat_list
  df_hnst <- df[1,]
  df_hest <- df[1,]
  df_name <- df[1,]
  df_pretty <- df[1,]
    
  # 1) Mean of the difference of means
  if (is.element("xbar_dm", stat_list)) {
    df$xbar_dm = rowMeans(x_b) - rowMeans(x_a)
    df_hnst$xbar_dm <- "<"
    df_hest$xbar_dm <- ">"
    df_pretty$xbar_dm <- "bar(x)[DM]"
  }
  
  # 2) Rel Means: mean of the difference in means divided by control group mean
  if (is.element("rxbar_dm", stat_list)) {
    df$rxbar_dm <- df$xbar_dm/rowMeans(x_a)
    df_hnst$rxbar_dm <- "<"
    df_hest$rxbar_dm <- ">"
    df_pretty$rxbar_dm <- "r*bar(x)[DM]"
  }
  
  # 3) Std of the difference in means
  if (is.element("sd_dm", stat_list)) {
    df$sd_dm <- sqrt(rowSds(x_a)^2/n_a + rowSds(x_b)^2/n_b)
    df_hnst$sd_dm <- "<"
    df_hest$sd_dm <- "<"
    df_pretty$sd_dm <- "s[DM]"
  }
  
  # 4) Relative STD: std of difference in means divided by midpoint between group means
  if (is.element("rsd_dm", stat_list)) {
    df$rsd_dm <- df$sd_dm / (rowMeans(x_a) + 0.5 * df$xbar_dm)
    df_hnst$rsd_dm <- "<"
    df_hest$rsd_dm <- "<"
    df_pretty$rsd_dm <- "r*s[DM]"
  }
  
  # 5) Bayes Factor
  if (is.element("bf", stat_list)) {
    df$bf <- row_bayesf_2s(x_a, x_b, paired = FALSE, deltas = deltas/df$sd_dm)
    # Delta is specified in units relative to standard deviation for BF
    df_hnst$bf <- "<"
    df_hest$bf <- ">"
    df_pretty$bf <- "BF"
  }
  
  # 6) NHST P-value: The more equal experiment will have a larger p-value
  # diff_z_score <- row_zscore_2s(x_b, x_a)
  if (is.element("pvalue", stat_list)) {
    df$pvalue <- row_ttest_2s(x_b, x_a, conf.level = conf.level)
    df_hnst$pvalue <- ">"
    df_hest$pvalue <- "<"
    df_pretty$pvalue <- "p[N]"
  }
  
  # 7) TOST p value (Two tailed equivalence test)
  if (is.element("tostp", stat_list)) {
    df$tostp <- row_tost_2s(x_b, x_a, deltas = deltas,
                            conf.level = conf.level)
    df_hnst$tostp <- "<"
    df_hest$tostp <- ">"
    df_pretty$tostp <- "p[E]"
  }
  
  # 8) Second Generation P value
  if (is.element("p2", stat_list)) {
    df$p2 <- row_sgpv(x_a, x_b, null.los = -deltas, null.his = deltas, conf.level = conf.level)
    df_hnst$p2 <- ">"
    df_hest$p2 <- "<"
    df_pretty$p2 <- "p[delta]"
  }
  
  # 8) Cohens D
  if (is.element("cohend", stat_list)) {
    df$cohend <- row_cohend(x_a, x_b)
    df_hnst$cohend <- "<"
    df_hest$cohend <- ">"
    df_pretty$cohend <- "Cd"
  }
  
  # 9) most difference in means
  if (is.element("mdm", stat_list)) {
    df$mdm <- row_mdm(x_a, x_b, conf.level = conf.level)
    df_hnst$mdm <- "<"
    df_hest$mdm <- ">"
    df_pretty$mdm <- "delta[M]"
  }
  
  # 10) Relative most difference in means
  if (is.element("rmdm", stat_list)) {
    df$rmdm <- row_rmdm(x_a, x_b, conf.level = conf.level) 
    df_hnst$rmdm <- "<"
    df_hest$rmdm <- ">"
    df_pretty$rmdm <- "r*delta[M]"
  }
  
  # 11) Least Difference in Means
  if (is.element("ldm", stat_list)) {
    df$ldm <- row_ldm(x_a, x_b, conf.level = conf.level)
    df_hnst$ldm <- "<"
    df_hest$ldm <- ">"
    df_pretty$ldm <- "delta[L]"
  }
  
  # 12) Relative Least Difference in Means
  if (is.element("rldm", stat_list)) {
    df$rldm <- row_rldm(x_a, x_b, conf.level = conf.level)
    df_hnst$rldm <- "<"
    df_hest$rldm <- ">"
    df_pretty$rldm <- "r*delta[L]"
  }
  
  # 13) Random group
  if (is.element("rand", stat_list)) {
    df$rand <- rowMeans(matrix(rnorm(n_samples * 50, mean = 0, sd = 1), 
                               nrow = n_samples, ncol = 50))
    df_hnst$rand <- "<"
    df_hest$rand <- ">"
    df_pretty$rand <- "Rnd"
  }
  
  
  # Store attributes within df
  attr(df,"user_attributes") <- c("user_attributes","hat", "varnames", "varnames_pretty")
  attr(df,"hnst") <- df_hnst
  attr(df,"hest") <- df_hest
  attr(df,"varnames") <- colnames(df)
  attr(df,"varnames_pretty") <- df_pretty

  return(df)
}

