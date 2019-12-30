

library(dplyr)
library(tibble)

gglabel_height <- function (var,line_num, pad_mult) {
  line_ypos <- pad_mult*(max(var)-min(var)) + max(var) + line_num * pad_mult*(max(var)-min(var))
}




rfnorm_swept_param <- function (pop_mu,pop_sigma,n_samples=1e4, n_obs=100) {
  # Given mu and sigma, one of them being a number and the other being a vector of values, generate a series of
  # simulations of a folded distribution, one for each of the elements in the vectored parameter. 
  # The other variable remains fixed.
  #  
  # Inputs
  #  mu: population mean for normal distribution, can be single number or vector
  #  sigma: population sd for normal distribution, can be a single number or vector
  #  n_samples: the number of samples drawn for each simulation
  #  n_obs: the number of observations per sample (N)
  #  n_sims: number of simulations, equal to the length of the vectored parameter.
  
  
  # Determine which parameter is vectored to sweep simulations over
  param_ind <-1
  if (length(pop_mu)<length(pop_sigma)) {param_ind <-2}
  param_names <- c("pop_mu","pop_sigma")
  
  # Ether sigma or mu are swept
  n_sims <- length(get(param_names[param_ind]))
  
  # Simulation results to be calculated
  # Rows: simulation #
  # Cols: swept mu value
  norm_samples_mean <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  norm_samples_sd <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  
  fnorm_samples_mean <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  fnorm_samples_sd <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  
  
  for (n in seq(1,n_sims,1)) {
    # Generate  samples for this simulation, calculate mean and sd for each, store in column
    obs_data <- matrix(rnorm(n_samples*n_obs, mean = pop_mu[min(n,length(pop_mu))], 
                             sd = pop_sigma[min(n,length(pop_sigma))]),
                       nrow = n_samples, byrow = TRUE)
    
    # Mean of samples for each trial
    norm_samples_mean[,n] <- t(apply(obs_data,1,mean))
    # Standard deviation of samples for each trial
    norm_samples_sd[,n] <- t(apply(obs_data,1,sd))
    
    # Mean of samples for each trial
    fnorm_samples_mean[,n] <- t(apply(abs(obs_data),1,mean))
    # Standard deviation of samples for each trial
    fnorm_samples_sd[,n] <- t(apply(abs(obs_data),1,sd))
  }
  
  
  # Collapse matrix of mu's and sd's into a data frame
  tbl_norm_mean <- as.table(norm_samples_mean)
  df_norm_mean <- as.data.frame(tbl_norm_mean) %>% 
    rename(index=Var1,!!param_names[param_ind] := Var2, "sample_mean" = Freq) %>%
    mutate(distr=as.factor('N'))
  
  tbl_norm_sd <-as.table(norm_samples_sd, keep.rownames = FALSE,row.names = NULL,responseName = y)
  df_norm_sd <- as.data.frame(tbl_norm_sd) %>% 
    rename(index=Var1,!!param_names[param_ind] := Var2, sample_sd = Freq) %>%
    mutate(distr=as.factor('N'))
  
  
  # Collapse matrix of mu's and sd's into a data frame
  tbl_fnorm_mean <- as.table(fnorm_samples_mean)
  df_fnorm_mean <- as.data.frame(tbl_fnorm_mean) %>% 
    rename(index=Var1,!!param_names[param_ind] := Var2, "sample_mean" = Freq) %>%
    mutate(distr=as.factor('F'))
  
  tbl_fnorm_sd <-as.table(fnorm_samples_sd, keep.rownames = FALSE,row.names = NULL,responseName = y)
  df_fnorm_sd <- as.data.frame(tbl_fnorm_sd) %>% 
    rename(index=Var1,!!param_names[param_ind] := Var2, sample_sd = Freq) %>%
    mutate(distr=as.factor('F'))
  
  
  # Merge mean and sd datatables, and convert population swept parameter to factor
  df_norm <- merge(df_norm_mean, df_norm_sd, by=c("index","distr", param_names[param_ind]))
  df_fnorm <- merge(df_fnorm_mean, df_fnorm_sd, by=c("index","distr", param_names[param_ind]))
  df_results <- rbind(df_norm, df_fnorm)
  df_results$param_names[param_ind] <- as.factor(param_names[param_ind])
  
}