

# p.adjust(p, method = p.adjust.methods, n = length(p))


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
library(tidyr)

# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
library(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
library(TOSTER)
library(docstring)

# Helper Functions
source("R/mmd.r")
RowVar <- function(x, ...) rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)

# Test if a value is within an interval (inclusive)
within_ci <- function(ci, x) x >= ci[1] & x <= ci[2]

error_test_codes <-function(is_error_rate_zero, is_error_rate_alpha) {
#' Assign codes for error rate whether null hypothesis is rejected
#' 
#' @description Test if proportion of mmds that are above mu (error rate) are equal to 0.05 
#' and 0.00, return code do delineate combinations of both results
#'  (0) E = neither,  (1) E == 0.00
#'  (2) E == 0.05,  (3) E == both
  assign_code <- function(z,a) { if (z & a) {value=3} else if(z & !a) {value = 2}
    else if(!z & a) {value = 1} else {value = 0}; return(value)
  }
  error_rate_codes <- matrix(mapply(assign_code, is_error_rate_zero, is_error_rate_alpha, SIMPLIFY = TRUE,
         USE.NAMES = FALSE), nrow = dim(is_error_rate_zero)[1], ncol = dim(is_error_rate_zero)[2])
  colnames(error_rate_codes) <- colnames(is_error_rate_zero)
  rownames(error_rate_codes) <- rownames(is_error_rate_zero)
  return(error_rate_codes)
}

stats_param_sweep <- function( mus, sigmas, n_samples, n_obs, out_path, 
                              mu_ov_sigmas = NULL, overrride = TRUE, rand.seed=0) {
  #' Perform parameter sweep with specified mus and sigmas
  #' 
  #' @description QUantifies stats of a series of simulations of normal random 
  #' samples (n_samples) with a specified number of observations (n_obs) with 
  #' the specified values for mu (mus) and sigma (sigmas) for a normal 
  #' distribution.
  #'  
  #'  @param mus vector of mu values
  #'  @param sigmas vector of sigma values
  #'  @param n_samples number of samples (collection of measurements to simulate
  #'   one experiment)
  #'  @param n_obs number of observations
  #'  @param out_path output path to store results of calculation so simulations
  #'   do not have to be rerun
  #'  @param mu_ov_sigmas 
  #'  
  #'  @return dataframe that stores input parameters 

  # Store results to disk since calculations are significant
  set.seed(rand.seed)
  if (!file.exists(out_path) || overrride) {
    n_mus = max(c(length(mus), length(mu_ov_sigmas)))
    # Matrix diff and error of mmd
    
    if (!is.null(mus)) {dimnames = list(sigmas,mus)
    } else {dimnames = list(sigmas,mu_ov_sigmas)}

    df <- tibble(
      sigma = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mu = matrix(0L, nrow = length(sigmas), ncol = n_mus),
      mu_over_sigma = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                             dimnames = dimnames),
      mean_diff_mmd_mu = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                dimnames = dimnames),
      mean_rdiff_mmd_mu = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                 dimnames = dimnames),
      mean_mmd_error = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                              dimnames = dimnames),
      error_rate_tests = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                dimnames = dimnames),
      p_val_mmd_eq_zero = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                 dimnames = dimnames),
      p_val_mmd_eq_alpha = matrix(0L, nrow = length(sigmas), ncol = n_mus, 
                                  dimnames = dimnames))
    # browser()
    # Generate random samples based on random mu values
    for (r in seq(1, length(sigmas), by = 1)) {
      # With a vector of mu/sigma and sigma known, calculate mu
      if (!is.null(mu_ov_sigmas)) { mus = mu_ov_sigmas * sigmas[r]}
      
      for (c in seq(1, length(mus), by = 1)) {
        
        df$sigma[r, c] = sigmas[r]
        df$mu[r, c] = mus[c]
        df$mu_over_sigma[r, c] = mus[c]/sigmas[r]
        
        # Get samples, where each row is a seperate sample, columns are observations
        x1 <- matrix(rnorm(n_samples * n_obs, mus[c], sigmas[r]), nrow = n_obs)
        # Calculate the mmd
        mmd = apply(x1, 2, function (x) mmd_normal_zdist(x, conf.level = 0.95) )
        
        # Difference mmd to mu
        df$mean_diff_mmd_mu[r, c] <- mean(mmd) - abs(mus[c])
        # Relative difference mmd to mu
        df$mean_rdiff_mmd_mu[r, c] <- df$mean_diff_mmd_mu[r, c] / sigmas[r]
        # Error rate mmd and mu
        n_successes <- sum(mmd < abs(mus[c]))
        df$mean_mmd_error[r, c] <- n_successes / n_samples
        # Caluate p-values for test against error rate == 0
        df$p_val_mmd_eq_zero[r, c] <- binom.test(
          n_successes, n_samples, p = 0.00, alternative = "two.sided", conf.level = 0.95)$p.value
        # Caluate p-values for test against error rate == alpha
        df$p_val_mmd_eq_alpha[r, c] <- binom.test(
          n_successes, n_samples, p = 0.05, alternative = "two.sided", conf.level = 0.95)$p.value
        # # Fraction of monte carlo tries where mmd was in error
        # df$mean_mmd_error[r, c] <- n_successes / n_samples
      }
    }
    # Replace INF with NaNs
    df$mean_rdiff_mmd_mu[!is.finite(df$mean_rdiff_mmd_mu)] <- NaN
    # Save an object to a file
    saveRDS(df, file = out_path)
  } else {
    # Restore the object
    df <- readRDS(file = out_path)
  }
  return(df)
}



row_locate_binary_bounds <- function (xb){
  
  # Helper function to flip matrix left to right
  fliplr <- function(x) x[,ncol(x):1]
  
  # Code assumes that TRUE is on the left portion of matrix, if FALSE is there 
  # instead then invert
  if (sum(xb[,1]) < dim(xb)[1]) {
    xb = !xb
  }
  
  
  x <- matrix(as.numeric(xb),ncol = ncol(xb), nrow = nrow(xb))
  a<-t(apply(x,1,cumsum))
  b<-fliplr(t(apply(fliplr(!x),1,cumsum)))
  
  
  c <- b+a
  
  
  row_max_inds <- apply(c,1, function(x) mean(which(x == max(x))))
  return(row_max_inds)
  
}


locate_bidir_binary_thresh <- function(mus = NULL, sigmas, n_obs,temp_name, mu_ov_sigmas = NULL, rand.seed=0){
  #' Locate all 4 error transition boundaries
  #' Positive and negative direction, 0 and alpha error
  set.seed(rand.seed)
  
  n_cols <- max(c(length(mus), length(mu_ov_sigmas)))
  if (!is.null(mus)) {
    col_values = mus; varname <- "critical_mu"; neg_mu_ov_sigmas = mu_ov_sigmas;
  } else {
    col_values = mu_ov_sigmas; varname <- "critical_mu_over_sigma"; neg_mu_ov_sigmas = -mu_ov_sigmas;
  }
  
  # Run simulations calculating error of mmd with mu and sigma swept
  df_right <- stats_param_sweep(mus = mus, sigmas = sigmas, 
                                n_samples, n_obs, paste("temp/right_", temp_name,".rds"),
                                mu_ov_sigmas) 
  df_right$side <- as.factor("Right")
  
  # Run simulations calculating error of mmd with mu and sigma swept
  df_left <- stats_param_sweep(mus = mus, sigmas = sigmas, n_samples, n_obs, 
                               paste("temp/right_", temp_name,".rds"), 
                               mu_ov_sigmas= neg_mu_ov_sigmas) 
  df_left$side <- as.factor("Left")
  
  
  # Equivalence test versus middle column of same row
  p_threshold = 0.05 /(n_cols*length(sigmas))
  df_zero <- rbind(
    tibble(er = "0", side = "right", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_right$p_val_mmd_eq_zero  < p_threshold)),
    tibble(er = "alpha", side = "right", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_right$p_val_mmd_eq_alpha < p_threshold)))
  # Convert index of error rate transition to mu/sigma value
  df_zero[[varname]] <- 
    approx(x=1:n_cols,y = abs(col_values),
           xout = df_zero$critical_val_ind, 
           n = length(n_cols)*2L-1L)$y
  
  df_alpha <- rbind(
    tibble(er = "0", side = "left", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_left$p_val_mmd_eq_zero  < p_threshold)),
    tibble(er = "alpha", side = "left", sigma = sigmas, 
           critical_val_ind = row_locate_binary_bounds(
             df_left$p_val_mmd_eq_alpha  < p_threshold)))
  # Convert index of error rate transition to mu/sigma value
  df_alpha[[varname]] <- 
    approx(x=1:n_cols,y = abs(col_values),
           xout = df_alpha$critical_val_ind, 
           n = length(n_cols)*2L-1L)$y
  
  # Concatenate right and left dataframes of results
  df_crit <- rbind(df_zero,df_alpha)
  
  df_crit[[varname]][df_crit$side=="left"] = -
    abs(df_crit[[varname]][df_crit$side=="left"])
  # * abs() added just in case we are running the code multiple times when debugging
  
  
  return(df_crit)
  
}



rowcol_pearson <- function(m, row_vals, col_vals) {
  #' Calculate pearson correlation on a heatmap row by row and column by column,
  #' return bonferroni corrected p values
  #' 
  bonf_corr = dim(m)[1] + dim(m)[2]
  
  row_pearson_p = rep(0, dim(m)[1])
  for (r in seq( dim(m)[1])) {
    row_pearson_p[r] = cor.test(unname(m[r,]),abs(col_vals),method = "pearson")$p.value/bonf_corr  
  }
  # print(row_pearson_p)
  col_pearson_p = rep(0, dim(m)[2])
  for (c in seq(dim(m)[2])){
    col_pearson_p[c] = cor.test(unname(m[,c]),row_vals,method = "pearson")$p.value/bonf_corr     
  } 
  # print(col_pearson_p)
  
  row_sig_labels <- ifelse(row_pearson_p<0.05, "*","")
  col_sig_labels <- ifelse(col_pearson_p<0.05, "*","")
  
  pear_stats <- list(row_sig_labels,col_sig_labels,row_pearson_p,col_pearson_p)
  names(pear_stats)<-c("row_sig_labels","col_sig_labels","row_pearson_p","col_pearson_p")
  return(pear_stats)
}





# General variables
fig_num = "3"
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)

n_samples <- 1e3
rand.seed <- 0





# 2D visualization of mmd difference and error rate over mu and sigma
#                                                                              #
#______________________________________________________________________________#
mus <- seq(-2.5, 2.5, by = .1)
sigmas <- seq(.1, 5, by = .1)
n_obs <- 50

# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(mus, sigmas, n_samples, n_obs, "temp/mmd_Error_2D_mu_vs_sigma.rds") 

# Assemble results into square matrix
df_pear <- rowcol_pearson(df_results$mean_diff_mmd_mu, sigmas, mus)
# Difference MMD from Mu
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_diff_mmd_mu)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  # annotate("text", x=max(mus), y=sigmas+.02, label = df_pear$row_sig_labels, size=2,vjust=1) +
  # annotate("text", x=mus, y=max(sigmas)+.02, label = df_pear$col_sig_labels, size=2,vjust=1) +
  scale_fill_gradientn(colors=c("blue","white", "#C00000"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(0,2)) +
  # coord_cartesian(xlim = c(min(mus)+0.1, max(mus)+0.1), # This focuses the x-axis on the range of interest
  #                 clip = 'off') +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a1 mmd diff.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 


# Difference MMD from Mu
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_rdiff_mmd_mu)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df_pear <- rowcol_pearson(df_results$mean_rdiff_mmd_mu, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  # annotate("text", x=max(mus), y=sigmas+.02, label = df_pear$row_sig_labels, size=2,vjust=1) +
  # annotate("text", x=mus, y=max(sigmas)+.02, label = df_pear$col_sig_labels, size=2,vjust=1) +
  scale_fill_gradientn(colors=c("blue","white", "#C00000"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal")) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a2 mmd rdiff.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 


# Difference MMD from Mu
# COnvert from matrix to dataframe
df <- cbind(sigma = sigmas, as_tibble(df_results$mean_mmd_error)) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df_pear <- rowcol_pearson(df_results$mean_mmd_error, sigmas, mus)
# Plot heatmap
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size=8) +
  # annotate("text", x=max(mus), y=sigmas+.02, label = df_pear$row_sig_labels, size=2,vjust=1) +
  # annotate("text", x=mus, y=max(sigmas)+.02, label = df_pear$col_sig_labels, size=2,vjust=1) +
  scale_fill_gradientn(colors=c("blue","white", "#C00000"), guide = guide_colorbar
                       (raster = T, frame.colour = c("black"), frame.linewidth = .5,
                         ticks.colour = "black",  direction = "horizontal"),
                       limits=c(0,.1)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.justification = "left",  legend.key.height = unit(.05, "inch"),
        legend.key.width = unit(.3, "inch"),legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(.1,"inch"))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_a3 mmd error rate 2D.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2.2,
          base_asp = 3, base_width = 2, dpi = 600) 




# Difference MMD from Mu
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$p_val_mmd_eq_zero > 0.05/(length(sigmas)*length(mu_ov_sigmas)), 
  df_results$p_val_mmd_eq_alpha > 0.05/(length(sigmas)*length(mu_ov_sigmas))))) %>% gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) + 
  geom_tile()+ 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "#C00000", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_b4 mmd error test.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600) 





# 2D visualization of mmd difference and error rate over sigma and mu/sigma
#                                                                              #
#______________________________________________________________________________#
sigmas <- seq(.1, 5, by = .1)
mu_ov_sigmas <- seq (-.5, .5, by=0.01)
n_obs <- 50
set.seed(rand.seed)
# Run simulations calculating error of mmd with mu and sigma swept
df_results <- stats_param_sweep(NULL, sigmas, n_samples, n_obs,
                           "temp/mmd_Error_2D_mu_over_sigma_vs_sigma.rds", mu_ov_sigmas)
# Difference MMD from Mu
# COnvert from matrix to data frame
df <- cbind(sigma = sigmas, as_tibble(error_test_codes(
  df_results$p_val_mmd_eq_zero > 0.05/(length(sigmas)*length(mu_ov_sigmas)),
  df_results$p_val_mmd_eq_alpha > 0.05/(length(sigmas)*length(mu_ov_sigmas))))) %>% 
  gather(mu, z, -sigma)
df$mu <- as.numeric(df$mu)
df$sigma <- as.numeric(df$sigma)
df$z <- factor(df$z,levels = c("0","1","2","3"))
# Plot heat map
gg<- ggplot(df, aes(mu, sigma, fill= z)) +
  geom_tile()+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression(mu/sigma)) + ylab(expression(sigma)) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=c("white", "#C00000", "blue","purple"),drop=FALSE) +
  geom_vline(xintercept=0, color="black", size=0.2) +
  theme(legend.position="none")
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_b5 mmd error test mu_over_sigma.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2, dpi = 600)






# Identify location of coverage error boundaries with mu space
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
sigmas <- seq(.1, 5, by = .1); 
mus <- seq (.1/5, 2, by=0.02)
mu_ov_sigmas = NULL

df_crit_mu <- locate_bidir_binary_thresh(mus = mus,sigmas = sigmas, n_obs = n_obs, 
                                      temp_name = "mmd_Error_mu_vs_sigma.rds",
                                      mu_ov_sigmas = mu_ov_sigmas, rand.seed = rand.seed)
df_crit_mu$merge = paste(df_crit_mu$er, df_crit_mu$side)

# Plot of each boundary separate
df_pearson <- df_crit_mu %>% group_by(er,side) %>% summarize(pearson = cor.test(
  critical_mu, sigma,method = "pearson")$p.value)
df_pearson$adj_pearson <- p.adjust(df_pearson$pearson,"bonferroni")
gg <- ggplot(data=df_crit_mu,aes(x=sigma, y=critical_mu,
                                 shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  xlab(expression(sigma)) + 
  ylab(expression(mu)) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))+
  # scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  ")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_d mmd boundaries over mu.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)

# Box and Whiskers of Coverage Error Transition Region
p <- ggplot(df_crit_mu, aes(x=er, color = side,#group = interaction(er, side), 
                         y = abs(critical_mu))) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9), outlier.shape = NA) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") + 
  ylab(expression(abs(~mu~phantom(.))))+
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
p
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_c mmd transition over mu.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(abs(critical_mu) ~ er + side, data = df_crit_mu)
summary_mu <- summary(res.aov2)
capture.output(summary_mu, file = paste("figure/F", fig_num, "/F", fig_num,
                                        "_c mmd transition over mu.txt", sep=""))


# Identify location of coverage error boundaries with mu/sigma space
#                                                                              #
#______________________________________________________________________________#
n_obs <- 50
mus <- NULL
sigmas <- seq(.1, 5, by = .1); 
mu_ov_sigmas <- seq (0.10, 0.40, by=0.001)

df_crit_mu_ov_sigma <- locate_bidir_binary_thresh(mus=NULL, sigmas = sigmas, n_obs = n_obs, 
                                      temp_name = "mmd_Error_mu_over_sigma_vs_sigma.rds", 
                                      mu_ov_sigmas = mu_ov_sigmas, rand.seed = rand.seed)
df_crit_mu_ov_sigma$merge = paste(df_crit_mu_ov_sigma$er, df_crit_mu_ov_sigma$side)

# Plot of each boundary separate
df_pearson <- df_crit_mu_ov_sigma %>% group_by(er,side) %>% summarize(pearson = cor.test(
  critical_mu_over_sigma, sigma,method = "pearson")$p.value)
df_pearson$adj_pearson <- p.adjust(df_pearson$pearson,"bonferroni")
gg <- ggplot(data=df_crit_mu_ov_sigma,aes(x=sigma, y=critical_mu_over_sigma,
                                          shape = er, color = side)) +
  # facet_grid(rows = vars(er)) +
  geom_point(size=1) +
  geom_hline(yintercept=0,size=0.5) +
  theme_classic(base_size = 8) +
  xlab(expression(sigma)) + 
  ylab(expression(mu*phantom(.)*phantom(.)/sigma)) + 
  theme(legend.position="none") +
  scale_color_manual(values = c("#66c2a5", "#fc8d62"))
  # scale_color_discrete(name = "", labels = c("-Ra  ","-R0 ","+R0  ","+Ra  "))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_d mmd boundaries over mu_sigma.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 2.5, dpi = 600)

# Box and Whiskers of Coverage Error Transition Region
p <- ggplot(df_crit_mu_ov_sigma, aes(x=er, color = side,#group = interaction(er, side), 
                         y = abs(critical_mu_over_sigma))) + 
  #geom_violin( position = position_dodge( width = 0.9)) + 
  geom_boxplot( width = 0.2,position = position_dodge( width = 0.9), outlier.shape = NA) +
  theme_classic(base_size = 8) + theme(legend.position="none", 
                                       axis.title.x = element_blank()) +
  xlab("Error Rate Null Hypothesis") +  
  ylab(expression(abs(~mu*phantom(.))*phantom(.)/sigma))+
scale_color_manual(values = c("#66c2a5", "#fc8d62"))
p
save_plot(paste("figure/F", fig_num, "/F", fig_num, "_c mmd transition over mu sigma.tiff", 
                sep = ""), p, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)

res.aov2 <- aov(abs(critical_mu_over_sigma) ~ er + side , data = df_crit_mu_ov_sigma)
summary_mu_ov_sigma <- summary(res.aov2)
capture.output(summary_mu_ov_sigma, file = paste("figure/F", fig_num, "/F", fig_num,
                                        "_c mmd transition over mu sigma.txt", sep=""))


