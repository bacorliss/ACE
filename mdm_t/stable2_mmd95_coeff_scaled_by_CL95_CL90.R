
#' Produce stat look-up tables for calculating MDM


# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(broom)
p_load(scales)
p_load(ggplot2)
p_load(dplyr)
p_load(plyr)
p_load(grid)
p_load(gridExtra)
p_load(colorspace)
p_load("RColorBrewer")
p_load(cowplot)
# User defined libraries
source("R/mdm.R")
source("R/row_stats_toolbox.R")

# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
out_dir = file.path(getwd(), paste(base_dir, "/figure/T2",sep=""))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mus = c(seq(0, 0.1,0.001), seq(0.11,0.5,0.01), 1, 2, 5, 10, 50, 100,1000)
sigmas = runif(length(mus),0.1, 2)
n_samples = 35
n_obs = 50
set.seed(0)
df_coeff <- data.frame(mu=mus, sigma=sigmas, mean_mdm_96 = rep(0,length(mus)),
                       sd_mdm_95 = rep(0,length(mus)), mean_mabs_cl_95 = rep(0,length(mus)),
                       sd_mabs_cl_95 = rep(0,length(mus)), mean_maabs_cl_90 = rep(0,length(mus)), 
                       sd_mabs_cl_90 = rep(0,length(mus)))
# Sample around mean
for (n in seq_along(mus)) {  # print(mus[n])
  # For each mu, generate samples, align them, calculate mean MDM, CI_95, CL_90
  xi <- matrix(rnorm(n_samples*n_obs,mean = mus[n],sd=sigmas), nrow = n_samples, byrow = TRUE)
  
  # Normalize samples (x_bar = mu and sd = 1)
  xnorm <- (xi - rowMeans(xi))/rowSds(xi) + mus[n]
  
  # Calculate MDM
  mdm_95 <- apply(xnorm, 1, mdm_tdist)
  df_coeff$mean_mdm_95[n] <-  mean(mdm_95)
  df_coeff$sd_mdm_95[n] <-    sd(mdm_95)
  # Calculate 90% max abs CL
  mabs_cl_90 <- apply(xnorm, 1, function (x)  max_abs_cl_mean_z(x=x, alpha=0.10) )
  df_coeff$mean_mabs_cl_90[n] <- mean(mabs_cl_90)
  df_coeff$sd_mabs_cl_90[n] <-   sd(mabs_cl_90)
  # Calculate 95% max abs CL
  mabs_cl_95 <- apply(xnorm, 1, function (x)  max_abs_cl_mean_z(x=x, alpha=0.05) )
  df_coeff$mean_mabs_cl_95[n] <- mean(mabs_cl_95)
  df_coeff$sd_mabs_cl_95[n] <-   sd(mabs_cl_95)
  
}
# # Calculate Coefficient for mdm
df_coeff$coeff_mdm_95 <- (df_coeff$mean_mdm_95-df_coeff$mean_mabs_cl_90) /
  (df_coeff$mean_mabs_cl_95 - df_coeff$mean_mabs_cl_90)


# Export LU table to disk
df_lut = data.frame(abs_nmu = df_coeff$mu, coeff_mdm_95 = df_coeff$coeff_mdm_95)
write.csv(x=df_lut, file=file.path(out_dir,"coeff_mdm_CLa_CL2a.csv"))
