# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)

p_load(broom)
p_load(scales)
p_load(ggplot2)
p_load(plyr)
p_load(dplyr)
p_load(grid)
p_load(gridExtra)
p_load(colorspace)
p_load("RColorBrewer")
p_load(cowplot)
p_load(boot)
p_load(tidyr)
p_load(rbenchmark)
# Calculate conifdence intervals
source("R/mmd.R")
source("R/row_effect_sizes.R")
source("R/norm_versus_fnorm.R")
base_dir = "mmd_z"

fig_num = "4"
dir.create(file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep="")), showWarnings = FALSE,recursive = TRUE)

# Colors for figure
# Orange: MMD, [255, 127, 0], #FF7F00
# CL90 Blue, [255, 151, 151], #FF9797
# CL95 Red, [147, 173, 221], #93ADDD
# rTOST95, Purple, [112, 48, 160], #7030A0
col_pal <- data.frame(mmd95 = "#FF7F00", CL95 ="#FFA7A7", CL90 = "#BFCFEB", rTOST95 = "#7030A0")

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

RootSpline1 <- function (x, y, y0 = 0, verbose = FALSE) {
  #' Given (x, y) data, find x where the linear interpolation crosses y = y0
  #' the default value y0 = 0 implies root finding
  #' since linear interpolation is just a linear spline interpolation
  #' the function is named RootSpline1
  if (is.unsorted(x)) { ind <- order(x); x <- x[ind]; y <- y[ind] }
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] <= 0)
  ## analytical root finding
  xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  ## make a plot?
  if (verbose) {plot(x, y, "l"); abline(h = y0, lty = 2); points(xr, rep.int(y0, length(xr)))
  }
  xr
}



# Compare MMD, CL95, CL90 through null region
#
#------------------------------------------------------------------------------
# Generate repeatable sample of random normal numbers
mu = 0
sigma = 1
sigma_range = 0.32
n_obs = 35
set.seed(0)
y <- rnorm(n_obs, mean = mu, sd = sigma)

# Mean Shift samples to zero
y <- y - mean(y)
# Define sweep for shift in sample mean from zero
df = tibble(x = seq(from = -sigma_range*sigma, to = sigma_range*sigma, by = 0.02),grp=1)
# Replicate samples so each of the x offsets can be applied to a seprate row
y_expanded <- matrix(y, nrow = length(df$x),ncol = length(y), byrow = TRUE)
y_expanded + dim(df$x[row(y_expanded)])
y_sweep <- sweep(y_expanded,1, df$x,'+')

# Calculate most hidden difference
df$mmd_95    <- apply(y_sweep, 1, mmd_normal_zdist)
# Max absolute confidence level
df$mcl_90   <- apply(y_sweep, 1, function (x)  
  max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), alpha=0.10) ) 
df$mcl_95  <- apply(y_sweep, 1, function (x)  max(abs( 
  max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), alpha=0.05) )))
# T test p value
df$ttest_p_val  <- apply(y_sweep, 1, function (x)  t.test(x)$p.value )
# df$mmd_95 - df$mcl_90 
x_critical <-  RootSpline1(x=df$x,y=df$ttest_p_val,y0 = 0.05)
ci_labels = c(bquote(max((~abs(CI[95])))~" "),
              bquote(max((~abs(CI[90])))), bquote(MMD[95]))

## Subplot A,B: MMD[a] transitions from CI[a] to CI[a/2] as the sample mean departs from zero 
g1A = ggplot(data=df,mapping = aes(x=x,y=mmd_95))
g1A <- g1A + theme_classic() +
  geom_rect(aes(xmin=-Inf, xmax=x_critical[1], ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_rect(aes(xmin=x_critical[2], xmax=Inf, ymin=-Inf, ymax=Inf), fill = lighten('black',0.9)) +
  geom_line(aes(x=x,y=mcl_90, col="CL_90"), size=.8, linetype="solid", color = col_pal$CL90) +
  geom_line(aes(x=x,y=mcl_95, col="CL_95"), size=.8, linetype="solid",color = col_pal$CL95) +
  geom_point(aes(col="MMD_95"), shape = 16,size = 1, color = col_pal$mmd95) +
  xlab("") + ylab("f(x)") +
  theme_classic(base_size=8) + theme(
    axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
    legend.position = "right",legend.justification = c(0, 0.5), legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)), plot.margin = unit(c(0,0,0,0),"mm"), legend.box = "vertical",
    legend.text=element_text(size=8),legend.key.size = unit(.5,"line"), legend.text.align=0,
    axis.title.y = element_text(size = 8)) 
  # scale_y_continuous(labels=function(x) sprintf("%.1f", x))
g1A 

# Plot p value of sample sweep
g1B_labels = c("Conf. Level  ","t-test","Crit. Region")
g1B = ggplot(data=df,mapping = aes(x=x,y=ttest_p_val))
g1B <- g1B +
  geom_rect(aes(xmin=-Inf, xmax = x_critical[1], ymin = -Inf, ymax = Inf, fill = lighten('black',0.9))) +
  geom_rect(aes(xmin = x_critical[2], xmax = Inf, ymin = -Inf, ymax = Inf, fill = lighten('black',0.9))) +
  geom_hline(aes(yintercept=0.05, col="Conf. Level"), color="red") + 
  geom_point(aes(col="t-test "), shape = 16, size = 0.5, color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression(bar(x))) + ylab("P-value") +
  # scale_y_continuous(name="p value", limits=c(-.1, 1.1),breaks=seq(0, 1, .5)) +
  # scale_shape_manual("Legend", labels=g1B_labels, values = c(2,2,NA)) +
  # scale_linetype_manual("Legend", labels = g1B_labels, values=c("dotted", "dotted",NA))+
  # scale_color_manual("Legend", labels = g1B_labels,values=c('red','black',"grey")) +
  scale_fill_manual("",values = lighten('black',0.9), guide = guide_legend(override.aes = list(alpha = 1))) + 
  theme_classic(base_size=8) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = c(0,.5,1)) + 
  coord_cartesian(ylim=c(0,1))
g1B
# Use cowplot to align subplots
cg = plot_grid(g1A, g1B, label_size = 12, ncol=1,rel_heights = c(.5,.5) )
cg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "2AB_MMD_vs_CI95.tiff",sep=""),
          cg, ncol = 1, nrow = 1, base_height = 1.25,
          base_asp = 3, base_width = 5, dpi = 600) # paper="letter"
graphics.off()




# MMD and rTOST compared to CL95 and CL90
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
             mean_mmd95 = rep(0, length(mus)), sd_mmd95 = rep(0, length(mus)), 
             mean_rtost95 = length(mus),   sd_rtost95 = rep(0, length(mus)),
             mean_coeff_mmd95 = rep(0, length(mus)), sd_coeff_mmd95 = rep(0, length(mus)), 
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
  # Calculate MMD across samples
  mmd95 <- apply(xr, 1, function(x) mmd_normal_zdist(x, conf.level = 0.95))
  df$mean_mmd95[n]  <- mean(mmd95)
  df$sd_mmd95[n]    <- sd(mmd95)
  # Coeff relative
  coeff_mmd95 <- (mmd95-macl90)/(macl95-macl90)
  df$mean_coeff_mmd95[n]  <- mean(coeff_mmd95)
  df$sd_coeff_mmd95[n]    <- sd(coeff_mmd95)
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

# Plot MMD and rTOST normlaized to CL95 and CL90
gg <- ggplot(data = df_plot,(aes(x=mu, y=mean_y))) +
  geom_hline(aes(yintercept = 0, linetype = "CL_90"), color = "#377eb8",  size = 1, alpha = 0.2) +
  geom_hline(aes(yintercept = 1, linetype = "CL_95"), color = "#e41a1c", size = 1, alpha = 0.2) +
  geom_point(aes(color = metric), size = 1) +
  scale_linetype_manual(name = "", labels = c(expression(Max(abs(CL[95]))),
                                              expression(Max(abs(CL[90])))),
                        values = c("solid", "solid")) +
  scale_color_manual(name="", labels = c( expression(MMD[95]), expression(rTOST[95])),
                     values=c("#ff7f00","#984ea3")) +
  xlab(expression(bar(x))) + ylab("Coeff. CL [90-95]") +
  theme_classic(base_size = 8) + theme(legend.position = "none")
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "F3a_RTOST_vs_MMD.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.25, base_width = 3.5, dpi = 600) 




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
# Bland altman of agreement between MMD algorithms
gg <- ggplot(df_compare, aes(x=means,y=rdiffs)) +
  geom_hline(yintercept = 1.96*sd(df_compare$rdiffs), color = "grey30", linetype="dashed", size=0.25) +
  geom_hline(yintercept = -1.96*sd(df_compare$rdiffs), color = "grey30", linetype="dashed", size=0.25) +
  geom_hline(yintercept = 0, color="grey30", size=0.25)+
  geom_point(size=0.05) +
  xlab("Mean") + 
  ylab("Rel. Diff.") +
  theme_classic(base_size=8)+
  scale_y_continuous(labels = function(x) scales::scientific(x,digits = 2)) 
  # geom_blank(aes(y = 1.1E-15)) +
  # geom_blank(aes(y = -1.1E-15))
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "g_BA rTOST95 vs MACL90.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.25,
          base_asp = 3, base_width = 2, dpi = 600) 







# MMD transition as LUT between CL95 and CL90 across post-normalized samples
#
# -----------------------------------------------------------------------------
# Generate 1000 samples, loop through different shifts, and quantify MMD, UCL_95, UCL_90
mus = c(seq(0-0.00001,0.33,0.00001), .5, 1, 2,5, 10, 20,50, 100, 500, 1000, 10000)
# mus = c(seq(0, 0.1,0.001), seq(0.11,0.5,0.01), 1, 2, 5, 10, 50, 100,1000)

sigmas = runif(length(mus),0.1, 2)
n_samples = 50
n_obs = 50
set.seed(0)
df_coeff <- data.frame(mu=mus, sigma=sigmas, mean_mmd_96 = rep(0,length(mus)),
                       sd_mmd_95 = rep(0,length(mus)), mean_mabs_cl_95 = rep(0,length(mus)),
                       sd_mabs_cl_95 = rep(0,length(mus)), mean_maabs_cl_90 = rep(0,length(mus)), 
                       sd_mabs_cl_90 = rep(0,length(mus)))
# Sample around mean
for (n in seq_along(mus)) {  # print(mus[n])
  # For each mu, generate samples, align them, calculate mean MMD, CI_95, CL_90
  xi <- matrix(rnorm(n_samples*n_obs,mean = mus[n],sd=sigmas), nrow = n_samples, byrow = TRUE)
  
  # Normalize samples (x_bar = mu and sd = 1)
  xnorm <- (xi - rowMeans(xi))/rowSds(xi) + mus[n]
  
  # Calculate MMD
  mmd_95 <- apply(xnorm, 1, mmd_normal_zdist)
  df_coeff$mean_mmd_95[n] <-  mean(mmd_95)
  df_coeff$sd_mmd_95[n] <-    sd(mmd_95)
  # Calculate 90% max abs CL
  mabs_cl_90 <- apply(xnorm, 1, function (x)  max_abs_cl_mean_z(x=x, alpha=0.10) )
  df_coeff$mean_mabs_cl_90[n] <- mean(mabs_cl_90)
  df_coeff$sd_mabs_cl_90[n] <-   sd(mabs_cl_90)
  # Calculate 95% max abs CL
  mabs_cl_95 <- apply(xnorm, 1, function (x)  max_abs_cl_mean_z(x=x, alpha=0.05) )
  df_coeff$mean_mabs_cl_95[n] <- mean(mabs_cl_95)
  df_coeff$sd_mabs_cl_95[n] <-   sd(mabs_cl_95)
  # Calcualte mmd coeff
  coeffs_mmd_95 <- (mmd_95 - mabs_cl_90)/ (mabs_cl_95 - mabs_cl_90)
  df_coeff$mean_coeff_mmd_95[n] <- mean(coeffs_mmd_95)
  df_coeff$sd_coeff_mmd_95[n] <- sd(coeffs_mmd_95)
  
}
# # Calculate Coefficient for mmd
# df_coeff$coeff_mmd_95 <- (df_coeff$mean_mmd_95-df_coeff$mean_mabs_cl_90) / 
#   (df_coeff$mean_mabs_cl_95 - df_coeff$mean_mabs_cl_90)
mean(df_coeff$sd_mmd_95)
mean(df_coeff$sd_mabs_cl_90)
mean(df_coeff$sd_mabs_cl_95)

# Plot look up table results
gg <- ggplot(data = subset(df_coeff,mu<.3),aes(x=mu, y=mean_coeff_mmd_95)) +
  geom_line(size=0.15) +
  geom_ribbon(aes(ymin = mean_coeff_mmd_95-1.96*sd_coeff_mmd_95,
                  ymax=mean_coeff_mmd_95+1.96*sd_coeff_mmd_95), fill = "grey") +
  xlab(expression(abs(phantom(.)*mu*phantom(.))*phantom(.)/sigma)) +
  ylab(expression(Coeff.~MMD[95])) +
  theme_classic(base_size=8)
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "2C_Coeff_mmd_CLa_CL2a.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.25,
          base_asp = 3, base_width = 2, dpi = 600) # paper="letter"
graphics.off()

# Export LU table to disk
df_lut = data.frame(abs_nmu = df_coeff$mu, coeff_mmd_95 = df_coeff$mean_coeff_mmd_95)
write.csv(x=df_lut, file=file.path(getwd(),"/R/coeff_mmd_CLa_CL2a.csv"))


# Test agreement with MMD lut to MMD root Bland altman
#
#-------------------------------------------------------------------------------
mus = seq(0,0.5,0.001)
n_samples = 1000
set.seed(0)
mus = runif(n_samples, -50,50)
sigmas = runif(n_samples,.1,10)
# n_obs = 50

# Sample around mean
x_samples = t(mapply(function(x,y) rnorm(n_obs, mean=x, sd=y),mus,sigmas, SIMPLIFY = TRUE))

# Load csv Look up table to convert to spline interp function
df_lut <- read.csv(file=file.path(getwd(),"/R/coeff_mmd_CLa_CL2a.csv"))
interp_fun = splinefun(x=df_lut$abs_nmu, y=df_lut$coeff_mmd_95, method="fmm",  ties = mean)

# Function to determine 95% MMD with LUT
mmd_95_lut <- function (x,interp_fun) {
  mabs_cl_90 <- max_abs_cl_mean_z(x=x, alpha=0.10)
  mabs_cl_95 <- max_abs_cl_mean_z(x=x, alpha=0.05)
  # Normalized mu
  abs_nmu = abs(mean(x)/sd(x))
  coeff_mmd <- interp_fun(abs_nmu)
  mmd_95 <- coeff_mmd * (mabs_cl_95 - mabs_cl_90) + mabs_cl_90
  
  return(mmd_95)
}
#
df_compare <- data.frame(mmd_root = apply(x_samples, 1, mmd_normal_zdist), 
                         mmd_lut  = apply(x_samples, 1, function(x) mmd_95_lut(x, interp_fun)))
df_compare$diffs = df_compare$mmd_root - df_compare$mmd_lut
df_compare$means = rowMeans(cbind(df_compare$mmd_root,df_compare$mmd_lut ))
# Bland altman of agreement between MMD algorithms
gg <- ggplot(df_compare, aes(x=means,y=diffs)) +
  geom_hline(yintercept = 1.96*sd(df_compare$diffs), color = "grey30", linetype="dashed", size=0.25) +
  geom_hline(yintercept = -1.96*sd(df_compare$diffs), color = "grey30", linetype="dashed", size=0.25) +
  geom_hline(yintercept = 0, color="grey30", size=0.25)+
  geom_point(size=0.05) +
  # xlab(expression((MMD[root]+MMD[lut])/2)) + 
  # ylab(expression(MMD[root]-MMD[lut])) +
  xlab("Mean") + 
  ylab("Diff.") +
  theme_classic(base_size=8)
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "g_BA MMD root vs MMD lut.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.25,
          base_asp = 3, base_width = 2, dpi = 600) 




# Speed benchmark between MMD algorithms
#
#-------------------------------------------------------------------------------
mmd_root_time = rep(0,100)
mmd_lut_time = rep(0,100)
for (n in 1:100) {
  results <- benchmark("mmd_root" = {
    mmd_root = apply(x_samples, 1, mmd_normal_zdist)
  },
  "mmd_lut" = {
    mmd_lut = apply(x_samples, 1, function(x) mmd_95_lut(x, interp_fun))
  },
  replications = 1,
  columns = c("user.self"))
  mmd_root_time[n] <- results$user.self[2]
  mmd_lut_time[n] <- results$user.self[1]
}
t.test(mmd_root_time, mmd_lut_time)
1- 60*mean(mmd_lut_time)/(60*mean(mmd_root_time))
df_speed <- tibble(x = as.factor(c(rep("MMD[root]",100),rep("MMD[lut]",100))),
                   mmd_root = c(mmd_root_time, mmd_lut_time)*60)
gg <- ggplot(data = df_speed,  aes(x=x, y=mmd_root)) + 
  geom_boxplot( outlier.size = 1) + theme_classic(base_size = 8) +
  ylab("Sec./1000 Runs") + xlab("Algorithm") +
  scale_x_discrete(limits = rev(levels(df_speed$x)),
                   labels = c('MMD[root]' = expression(MMD[root]),
                              'MMD[lut]'   = expression(MMD[lut])))
gg
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "g_Speed MMD root vs MMD lut.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.25,
          base_asp = 3, base_width = 2, dpi = 600) 




## Explore trends between MMD and CI at a population level
#
#--------------------------------------------------------------------------------------#

# Generate 1000 samples, loop through different shifts, and quantify MMD, UCL_95, UCL_90
mus = seq(-1,1,.2)
sigma = 1
n_samples = 1000
n_obs = 35
set.seed(1)
# Sample around mean


# Make list of data frames to be concatenated at end of computations
df_list <- list()
for (n in seq(1,length(mus),1)) {
  
  y_sweep <- matrix(rnorm(n_samples*n_obs,mus[n],sigma), nrow = n_samples, byrow = TRUE)
  # Calculate MMD and max abs confidence intervals
  mmd_95    <- apply(y_sweep, 1, mmd_normal_zdist)
  # Two tailed confidence intervals
  mcl_90   <- apply(y_sweep, 1, function (x)  
    max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), alpha=0.10) )
  mcl_95  <- apply(y_sweep, 1, function (x)  
    max_abs_cl_mean_z(mean(x), sd(x)/sqrt(length(x)), alpha=0.05) )
  # mcl_90_t   <- apply(y_sweep, 1, function (x)  max(abs(t.test(x, conf.level = 1-0.10)$conf.int )))
  # mcl_95_t  <- apply(y_sweep, 1, function (x)  max(abs(t.test(x, conf.level = 1-0.05)$conf.int )))
  mmd_diff <- mmd_95 - mcl_90
  ci_diff <- mcl_95 - mcl_90
  coeff_mmd95 <- mmd_diff/ci_diff
  print(mus[n])
  df_list[[n]] = tibble(n = as.factor(n), 
                        coeff_mmd95 = coeff_mmd95, mmd_95 = mmd_95, 
                        mcl_90 = mcl_90, mcl_95 = mcl_95)
  df_list[[n]]$mu <-as.factor(mus[n])
}

df <- ldply(df_list, rbind)
# df$mu <- as.factor(df$mu)
# Get groups means and CI of mean
df_plotted <- df %>% group_by(mu) %>% 
  summarize(mean_coeff_mmd95 = mean(coeff_mmd95)) 
            # lcl_mu = boot.ci(boot(df$coeff_mmd95, statistic=
            #                         function(data, i) {return(mean(data[i]))}, R=1000), 
            #                  conf=0.95, type="bca")$bca[4],
            # ucl_mu = boot.ci(boot(df$coeff_mmd95, statistic=
            #                         function(data, i) {return(mean(data[i]))}, R=1000), 
            #                  conf=0.95, type="bca")$bca[5])
df_plotted$mu <- as.factor(df_plotted)
# Pairwise test between groups
pw_means <- pairwise.t.test(df$coeff_mmd95, df$mu, p.adjust.method = "bonferroni",
                                paired = TRUE, alternative = "two.sided")
adj_sig_str <- adjacent_compare_str(pw_means$p.value<0.05,'*')
  
# Plotting MMD vs coefficient population level
gg <- ggplot(df, aes(x = mu, y=coeff_mmd95)) +
  geom_violin(scale = "width", fill = "grey85", color="white") +
  geom_point(data=df_plotted, aes(x=mu,y=mean_coeff_mmd95), color ="#FF7F00") +
  geom_hline(aes(yintercept = 0), color = "#377eb8",  size = 1, alpha = 0.2) +
  geom_hline(aes(yintercept = 1), color = "#e41a1c", size = 1, alpha = 0.2) +
  xlab(expression(mu)) + ylab("Coeff. CL [90, 95]") +
  geom_text(data = data.frame(), aes(x = as.factor(mus), y = rep(1.03, length(mus)), 
                label=adj_sig_str),  size=5, position = position_nudge(x = 0.5, y = 0), color="black")+
  coord_cartesian(ylim=c(0,1.05))+
  theme_classic(base_size=8)
gg
# Export figure to disk
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "5a_coeff_MMD_from_pop.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 5, dpi = 600) # paper="letter"
graphics.off()



