


#' Repeated samples from swept normal population parameters showing the the rev_tost
#' is equivalent to the max( abs( 90% CI))

# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# Load packages
p_load(ggplot2)
p_load(tibble)
p_load(equivalence)
p_load(TOSTER)
p_load(dplyr)
p_load(tidyr)
p_load(cowplot)
# User defined functions
source("R/row_stats_toolbox.R")
source("R/mdm.R")



# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
fig_num = "3"
dir.create(file.path(getwd(), paste(base_dir,"/figure/SF",fig_num,sep="")), 
           showWarnings = FALSE, recursive = TRUE)
rand.seed = 0;

ztest <- function(x, y = NULL, mu = 0, conf.level= 0.95, alternative = "two.sided") {
  #' Classic Z test of one or two samples
  #' 
  #' 
  #' 
  #' 
  if (is.null(y)) { # One sample
    n_x <- length(x)
    z_stat <- (mean(x) - mu) / (sd(x)/sqrt(n_x))
  } else {          # Two sample
    n_x <- length(x); n_y = length(y)
    z_stat <- (mean(y) - mean(x)- mu) / sqrt(sd(x)^2/n_x - sd(y)^2/n_y)
  }
  
  two_tail_p   <- 2* pnorm(abs(z_stat),lower.tail = FALSE, log.p = FALSE)
  lesser_p <- pnorm(z_stat,lower.tail = TRUE, log.p = FALSE)
  greater_p <- 1 - pnorm(z_stat,lower.tail = TRUE, log.p = FALSE)
  
  if (alternative=="two.sided") {# Two_tail (not equal to): H0 u = u0, H1: u != u0
    p <- two_tail_p
  } else if (alternative=="lesser") {# Lower Tail (less than): H0 u >= u0, H1: u < u0
    p <- lesser_p
  } else if (alternative=="greater") {# Upper Tail (greater than): H0 u <= u0, H1: u > u0
    p <- greater_p 
  } else {errorCondition("ztest: unknown alternative hypothesis")}
  
  # df <- tibble(two_tail_p=two_tail_p, lesser_p = lesser_p, greater_p = greater_p)
  return(p)
}

ztost <- function(x, delta) {
  #' Two one sided Tests for z distribution, used for equivalence testing
  #' 
  #' 
  # Equivalence: -delta <=  u1 < delta
  # z_low  <- (mean(x) - delta)/(sd(x) / sqrt(length(x)))
  # z_high <- (mean(x) + delta)/(sd(x) / sqrt(length(x)))

  p_lower <- ztest(x, mu = -abs(delta), alternative = "greater")
  p_upper <- ztest(x, mu =  abs(delta), alternative = "lesser" )
  # H01: delta <= -delta_0,     Ha1: delta >  -delta_0
  # p1 <- pnorm(-zp_lower_low)
  # H02: delta >= delta_0,      Ha2: delta <  delta_0
  # p2 <- pnorm(z_high)
  
  p = max(c(p_lower,p_upper))
  return(p)
}

rztost <- function(xs, alpha) {
  #' Calculates the reverse two one sided tests for z distribution, returns the 
  #' largest equivalence region where results are still significant (p=0.05)
  rtost_crit <- uniroot(function(x) ztost(xs, delta=x) - alpha,
          interval = c(0, abs(mean(xs)) + 6*sd(xs) ), 
          tol = .Machine$double.eps)$root
  return(rtost_crit)
}



# MDM and rTOST compared to CL95 and CL90
#
#-------------------------------------------------------------------------------
# n_sims = 101
n_samples = 1E3
n_obs = 50
mus <-  seq(-0.26, 0.26, by = .02)
sigmas <- rep(.1, length(mus))
# Initialize dataframe of metrics
df <- tibble(mu = mus, sigma = sigmas, 
             mean_macl90 = length(mus),          sd_macl90 = rep(0, length(mus)),
             mean_macl95 = length(mus),          sd_macl95 = rep(0, length(mus)),
             mean_mdm95 = rep(0, length(mus)), sd_mdm95 = rep(0, length(mus)), 
             mean_rtost95 = length(mus),   sd_rtost95 = rep(0, length(mus)),
             mean_coeff_mdm95 = rep(0, length(mus)), sd_coeff_mdm95 = rep(0, length(mus)), 
             mean_coeff_rtost95 = length(mus),   sd_coeff_rtost95 = rep(0, length(mus)))

x0 = t(sapply(1:n_samples, function(x) rnorm(n_obs, mean = 0, sd = 1),
                         simplify = TRUE))
x0 = (x0-rowMeans(x0))/rowSds(x0)

# For each simulation, draw samples from specified population parameters
for (n in seq_along(mus)) {
  xr = x0 + mus[n]
  # Calculate max( |CL_90| )
  macl90 <- apply(xr, 1, function(x) max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), a = 0.10))
  df$mean_macl90[n]  <- mean(macl90)
  df$sd_macl90[n]    <- sd(macl90)
  # Calculate max( |CL_95| )
  macl95 <- apply(xr, 1, function(x) max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), a = 0.05))
  df$mean_macl95[n]  <- mean(macl95)
  df$sd_macl95[n]    <- sd(macl95)  
  
  # Calculate MDM across samples
  mdm95 <- apply(xr, 1, function(x) mdm_tdist(x, conf.level = 0.95))
  df$mean_mdm95[n]  <- mean(mdm95)
  df$sd_mdm95[n]    <- sd(mdm95)
  # Coeff relative
  coeff_mdm95 <- (mdm95-macl90)/(macl95-macl90)
  df$mean_coeff_mdm95[n]  <- mean(coeff_mdm95)
  df$sd_coeff_mdm95[n]    <- sd(coeff_mdm95)
  
  # calculate reverse TOST across all samples
  rtost95 <- apply(xr, 1, function(x) rztost(x, alpha = 0.05)  )
  df$mean_rtost95[n]  <- mean(rtost95)
  df$sd_rtost95[n]    <- sd(rtost95)
  # Coeff relative
  coeff_rtost95 <- (rtost95-macl90)/(macl95-macl90)
  df$mean_coeff_rtost95[n] <- mean(coeff_rtost95)
  df$sd_coeff_rtost95[n]   <- sd(coeff_rtost95)
}

df_plot <- df %>% gather(metric, mean_y, starts_with("mean_coeff")) #%>%
  # gather(metric, sd_y, starts_with("sd_coeff"))
df_plot$metric <- as.factor(df_plot$metric)

# Plot MDM and rTOST normlaized to CL95 and CL90
gg <- ggplot(data = df_plot,(aes(x=mu, y=mean_y))) +
  geom_hline(aes(yintercept = 0, linetype = "CL_90"), color = "#e41a1c",  size = 2, alpha = 0.2) +
  geom_hline(aes(yintercept = 1, linetype = "CL_95"), color = "#377eb8", size = 2, alpha = 0.2) +
  geom_point(aes(color = metric), size = 1) +
  scale_linetype_manual(name = "", labels = c(expression(Max(abs(CL[95]))),
                                            expression(Max(abs(CL[90])))),
  values = c("solid", "solid")) +
  scale_color_manual(name="", labels = c( expression(MDM[95]), expression(rTOST[95])),
                     values=c("#ff7f00","#984ea3")) +
  xlab(expression(bar(x))) + ylab("Coeff. CL [90-95]") +
  theme_classic(base_size = 8)
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "F4_RTOST_vs_MDM.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.75, base_width = 4, dpi = 600) 



# Agreement between rTOST and MACL90
#
#-------------------------------------------------------------------------------
set.seed(rand.seed)
n_samples = 1E2
mus <-  seq(-2, 2, by = 0.1)
sigmas <- runif(n_samples, 0.1, 1)

df_list <- list() 
for (n in seq_along(mus)) {
  xr <- t(sapply(sigmas, function(x) rnorm(n_obs, mean = mus[n], sd = x),
                              simplify = TRUE))
  dfn <- tibble(mu = mus[n],sigma = sigmas, 
                macl90 = apply(xr, 1, function(x) max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), a = 0.10)),
                rtost95 = apply(xr, 1, function(x) rztost(x, alpha = 0.05)  ))
  df_list[[n]] <- dfn
}
df <- head(do.call(rbind, df_list))
df <- bind_rows(df_list, .id = "column_label")

df_compare <- tibble(mu = df$mu, sigma = df$sigma, 
                     rdiffs = (df$macl90 - df$rtost95)/rowMeans(cbind(df$macl90,df$rtost95 )),
                     means = rowMeans(cbind(df$macl90,df$rtost95 )))
# Bland altman of agreement between MDM algorithms
gg <- ggplot(df_compare, aes(x=means,y=rdiffs)) +
  geom_hline(yintercept = 1.96*sd(df_compare$rdiffs), color = "red", linetype="dashed", size=0.25) +
  geom_hline(yintercept = -1.96*sd(df_compare$rdiffs), color = "red", linetype="dashed", size=0.25) +
  geom_hline(yintercept = 0, color="blue", size=0.25)+
  geom_point(size=0.05) +
  xlab("Mean") + 
  ylab("Rel. Diff.") +
  theme_classic(base_size=8)+
  # scale_y_continuous(labels = scales::number_format(accuracy = 1e-15))
  geom_blank(aes(y = 1.1E-15))+
  geom_blank(aes(y = -1.1E-15))
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "g_BA rTOST95 vs MACL90.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 2, dpi = 600) 




