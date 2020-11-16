
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
source("R/mmd.R")
source("R/row_effect_sizes.R")



dir.create(file.path(getwd(), paste("figure/T2",sep="")), showWarnings = FALSE)

mus = c(seq(0, 0.1,0.001), seq(0.11,0.5,0.01), 1, 2, 5, 10, 50, 100,1000)

sigmas = runif(length(mus),0.1, 2)
n_samples = 35
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
  
}
# # Calculate Coefficient for mmd
df_coeff$coeff_mmd_95 <- (df_coeff$mean_mmd_95-df_coeff$mean_mabs_cl_90) /
  (df_coeff$mean_mabs_cl_95 - df_coeff$mean_mabs_cl_90)


# Export LU table to disk
df_lut = data.frame(abs_nmu = df_coeff$mu, coeff_mmd_95 = df_coeff$coeff_mmd_95)
write.csv(x=df_lut, file=file.path(getwd(),"figure/T2/coeff_mmd_CLa_CL2a.csv"))
