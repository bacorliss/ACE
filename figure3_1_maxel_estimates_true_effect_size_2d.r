



# Load packages
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)
library(rlang)
library(colorspace)
library(VGAM)
library(boot)
library(dplyr)
library(cowplot)
library(binom)
library(VGAM)
library(gplots)
library(RColorBrewer)
# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
library(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
library(TOSTER)

source("R/mhd.r")

# Helper Functions
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
within_ci <- function(ci,x) x>=ci[1] & x <=ci[2]



fig_basename = "f_3"
n_samples <- 1e3
rand_seed <-0


# 2D visualization of MHD differecne and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
# n_samples <- 1e4
n_obs <- 50

# Matrix diff and error of MHD
df_2d <- tibble(mean_diff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = length(mus)),
                mean_rdiff_mhd_mu = matrix(0L, nrow = length(sigmas), ncol = length(mus)),
                mean_mhd_error = matrix(0L, nrow = length(sigmas), ncol = length(mus)),
                error_rate_tests = matrix(0L, nrow = length(sigmas), ncol = length(mus)))



# Store results to disk since calculations are significant
set.seed(rand_seed)
if (!file.exists("temp/MHD_Error_2D.rds")) {
  
  # Generate random samples based on random mu values
  for (c in seq(1, length(mus), by = 1)) {
    for (r in seq(1, length(sigmas), by = 1)) {
      # Get samples, where each row is a seperate sample, columns are observations
      x1 <-
        matrix(rnorm(n_samples * n_obs, mus[c], sigmas[r]), nrow = n_obs)
      
      # Calculate the MHD
      mhd = apply(x1, 2, function (x) mhd_1sample_normal(x, alpha = 0.05))
      
      # Difference MHD to mu
      df_2d$mean_diff_mhd_mu[r, c] <- mean(mhd) - abs(mus[c])
      # Relative difference MHD to mu
      df_2d$mean_rdiff_mhd_mu[r, c] <- df_2d$mean_diff_mhd_mu[r, c] / sigmas[r]
      # Error rate MHD and mu
      n_successes <- sum(mhd < abs(mus[c]))
      df_2d$mean_mhd_error[r, c] <- n_successes / n_samples
      
      # Test if binomial coeefficient of error rate is equal to 1 or alpha
      error_rate_ci <-
        binom.test(sum(mhd < abs(mus[c])), n_samples, p = 0.05, 
                   alternative = "two.sided", conf.level = 0.95)[4][[1]]
      is_error_rate_p05 <- within_ci(error_rate_ci, 0.05)
      is_error_rate_0 <- within_ci(error_rate_ci,   0)
      if (is_error_rate_p05 & is_error_rate_0) {
        df_2d$error_rate_tests[r, c] = 3
      } else if (is_error_rate_0) {
        df_2d$error_rate_tests[r, c] = 2
      } else if (is_error_rate_p05) {
        df_2d$error_rate_tests[r, c] = 1
      }

      df_2d$mean_mhd_error[r, c] <- sum(mhd < abs(mus[c])) / n_samples
    }
  }
  
  # Replace INF with NaNs
  df_2d$mean_rdiff_mhd_mu[!is.finite(df_2d$mean_rdiff_mhd_mu)] <- NaN
  # Save an object to a file
  saveRDS(df_2d, file = "temp/MHD_Error_2D.rds")
  
  
} else {
  # Restore the object
  df_2d <- readRDS(file = "temp/MHD_Error_2D.rds")
}



# Visualize difference between MHD and mu
# https://sebastianraschka.com/Articles/heatmaps_in_r.html
png(paste("figure/", fig_basename, "a1 MHD diff.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white","blue", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, .4, length=100), seq(0.41, 0.9, length=100),seq(0.91,2,length=100))             

heatmap.2(df_2d$mean_diff_mhd_mu, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow=rep("",length(sigmas)), labCol= round(mus,1),
          col = my_palette,
          density.info='none', scale="none",
          cexRow = 1, cexCol = 1,
          denscol="black", keysize=1,
          key.par=list(mar=c(3.5,0,3,0)),
          breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          key.title = expression(paste("Difference",~MHD~-abs(phantom(.)*mu*phantom(.)))),
          na.color = "black",
          key.xlab = "",
          main = NULL,
          xlab(expression(mu)),
          ylab(expression(sigma)))
dev.off()



# Visualize relative difference between MHD and mu
png(paste("figure/", fig_basename, "a2 MHD rdiff.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "red"))(n = 199)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, .3, length=100), seq(0.31, .4, length=100))             

heatmap.2(df_2d$mean_rdiff_mhd_mu, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow = rep("",length(sigmas)), labCol= round(mus,1),
          col = my_palette,
          density.info = 'none', scale="none",
          cexRow = 1, cexCol = 1,
          denscol = "black", keysize=1,
          key.par = list(mar=c(3.5,0,3,0)),
          breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          key.title = expression(paste("Relative Difference",~(MHD~-abs(~mu))/~sigma)),
          na.color = "black",
          key.xlab = "",
          main = NULL,
          xlab(expression(mu)),
          ylab(expression(sigma)))
dev.off()



# Visualize error rate of MHD
png(paste("figure/", fig_basename, "a3 MHD error rate 2D.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, 0.039, length=100), seq(0.04, 0.06, length=100),seq(0.061, 0.08,length=100)) 

color_cull <- function(x) x[seq(1,round(length(x)*.73), by = 1)]
heatmap.2(df_2d$mean_mhd_error, Rowv = FALSE, Colv = FALSE, trace = "none", dendrogram = "none", 
          labRow = rev(round(sigmas,1)), labCol = round(mus,1),
          col = my_palette,
          density.info='none', scale="none",
          cexRow = 1, cexCol = 1, denscol="black", keysize=1, 
          key.par = list( mar = c(3.5,0,3,0)),
          lmat = rbind(c(5, 4, 2), c(6, 1, 3)), margins = c(3,0),
          lhei = c(2, 6), lwid = c(1, 10, 1),
          breaks = col_breaks,
          key.title = "",
          key.xlab = "",
          main = NULL, 
          xlab( expression(mu)), 
          ylab( expression(sigma)))
dev.off()

# 2D Heatmap of error (1) p == 0.05), (2) p == 1, (3) p == 0.05 & p == 1
png(paste("figure/", fig_basename, "a4 MHD error test.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white","blue", "red", "purple"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(0, .5, length=100), 
#                seq(0.51, 1.5, length=100),
#                seq(1.51,2,length=100))             

heatmap.2(df_2d$error_rate_tests, Rowv=FALSE, Colv=FALSE, trace="none", dendrogram = "none", 
          labRow=rev(round(sigmas,1)), labCol= round(mus,1),
          col = my_palette,
          density.info='none', scale="none",
          cexRow = 1, cexCol = 1,
          denscol="black", keysize=1,
          key = FALSE,
          # breaks = col_breaks,
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), margins=c(3,0),
          lhei=c(2, 6), lwid=c(1, 10, 1),
          # key.title = expression(paste("Difference",~MHD~-abs(phantom(.)*mu*phantom(.)))),
          na.color = "black",
          # key.xlab = "",
          main = NULL,
          xlab(expression(mu)),
          ylab(expression(sigma)))
dev.off()








# 2D visualization of MHD difference and error rate over sigma and mu/sigma
#                                                                              #
#______________________________________________________________________________#
# mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
mu_ov_sigmas <- seq (-.5, .5, by=0.01)
# n_samples <- 1e4
n_obs <- 50

# Matrix diff and error of MHD
df_2d <- tibble(mean_diff_mhd_mu = matrix(0L, nrow = length(sigmas), 
                                          ncol = length(mu_ov_sigmas)),
                mean_rdiff_mhd_mu = matrix(0L, nrow = length(sigmas), 
                                           ncol = length(mu_ov_sigmas)),
                mean_mhd_error = matrix(0L, nrow = length(sigmas), 
                                        ncol = length(mu_ov_sigmas)),
                error_rate_tests = matrix(0L, nrow = length(sigmas), 
                                          ncol = length(mu_ov_sigmas)))


# Store results to disk since calculations are significant
set.seed(rand_seed)
if (!file.exists("temp/MHD_Error_2D_mu_ov_sigmas.rds")) {
  
  # Generate random samples based on random mu values
  for (r in seq(1, length(sigmas), by = 1)) {
    # With a vector of mu/sigma and sigma known, calculate mu
    mus = mu_ov_sigmas * sigmas[r]
    
     for (c in seq(1, length(mu_ov_sigmas), by = 1)) {
    
      
      
      # Get samples, where each row is a seperate sample, columns are observations
      x1 <-
        matrix(rnorm(n_samples * n_obs, mus[c], sigmas[r]), nrow = n_obs)
      
      # Calculate the MHD
      mhd = apply(x1, 2, function (x)
        mhd_1sample_normal(x, alpha = 0.05))
      
      # Difference MHD to mu
      df_2d$mean_diff_mhd_mu[r, c] <- mean(mhd) - abs(mus[c])
      # Relative difference MHD to mu
      df_2d$mean_rdiff_mhd_mu[r, c] <- df_2d$mean_diff_mhd_mu[r, c] / sigmas[r]
      # Error rate MHD and mu
      n_successes <- sum(mhd < abs(mus[c]))
      df_2d$mean_mhd_error[r, c] <- n_successes / n_samples
      
    }
  }
  
  # Replace INF with NaNs
  df_2d$mean_rdiff_mhd_mu[!is.finite(df_2d$mean_rdiff_mhd_mu)] <- NaN
  # Save an object to a file
  saveRDS(df_2d, file = "temp/MHD_Error_2D_mu_ov_sigmas.rds")
  
  
} else {
  # Restore the object
  df_2d <- readRDS(file = "temp/MHD_Error_2D_mu_ov_sigmas.rds")
}







# Visualize error rate of MHD
png(paste("figure/", fig_basename, "a5 MHD error rate 2D.png"),    
    width = 2.3*300, height = 3*300, res = 300, pointsize = 8)       
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, 0.039, length=100), seq(0.04, 0.06, length=100),seq(0.061, 0.08,length=100)) 

color_cull <- function(x) x[seq(1,round(length(x)*.73), by = 1)]
heatmap.2(df_2d$mean_mhd_error, Rowv = FALSE, Colv = FALSE, trace = "none", dendrogram = "none", 
          labRow = rev(round(sigmas,1)), labCol = round(mu_ov_sigmas,2),
          col = my_palette,
          density.info='none', scale="none",
          cexRow = 1, cexCol = 1, denscol="black", keysize=1, 
          key.par = list( mar = c(3.5,0,3,0)),
          lmat = rbind(c(5, 4, 2), c(6, 1, 3)), margins = c(3,0),
          lhei = c(2, 6), lwid = c(1, 10, 1),
          breaks = col_breaks,
          key.title = "",
          key.xlab = "",
          main = NULL, 
          xlab( expression(mu)), 
          ylab( expression(sigma)))
dev.off()


# Grab mean error rates for mu==0
RowVar(t(df_2d$mean_mhd_error))
rowMeans(t(df_2d$mean_mhd_error))







# For each row
for (r in seq(1, length(sigmas), by = 1)) {
  mhd_error_rate_zero <- df_2d$mean_mhd_error[r,(length(mu_ov_sigmas)+1)/2]
  mhd_error_rate_end <- df_2d$mean_mhd_error[r,length(mu_ov_sigmas)]
  
  
  for (c in seq(1, length(mu_ov_sigmas), by = 1)) {

    # a <- TOSTtwo.raw(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,low_eqbound=-0.05/20, 
    #                  high_eqbound=0.05/20, alpha = 0.05, var.equal=FALSE,plot = FALSE, verbose = FALSE)
    
  }
}
# Equivalence test versus middle column of same row

# Equivalence test versus last column of same row





plot(rowMeans(t(df_2d$mean_mhd_error)))

ttest(df_2d$mean_mhd_error)

