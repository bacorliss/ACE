
# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(dplyr)
p_load(tibble)
p_load(broom)

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
  param_names <- c("pop_mu","pop_sigma")
  param_ind <-1
  if (length(pop_mu)<length(pop_sigma)) {param_ind <-2}

  # Ether sigma or mu are swept
  n_sims <- length(get(param_names[param_ind]))
  
  # Simulation results to be calculated
  # Rows: simulation #
  # Cols: swept mu value
  # browser();
  norm_samples_mean <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  norm_samples_sd <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  
  fnorm_samples_mean <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  fnorm_samples_sd <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
  
  
  for (n in seq(1,n_sims,1)) {
    # Generate  samples for this simulation, calculate mean and sd for each, store in column
    obs_data <- matrix(rnorm(n_samples*n_obs, mean = pop_mus[min(n,length(pop_mus))], 
                             sd = pop_sigmas[min(n,length(pop_sigmas))]),
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
  
  # Lookup table to convert levels of swept variable from chars to actual numeric values
  lu_table <- get(param_names[param_ind])
  df_results[[(param_names[param_ind])]] <- as.factor(unname(lu_table[as.numeric(df_results[[(param_names[param_ind])]])]))
  
  # browser();
  return(df_results)
}


norm_fnorm_stats <- function (df_results, ind_varname) {

  
  # ANOVA for N
  summary(aov(sample_mean ~ get(ind_varname), data = subset(df_results, distr=='N')))
  summary(aov(sample_sd ~ get(ind_varname), data = subset(df_results, distr=='N')))
  # ANOVE for FN
  summary(aov(sample_mean ~ get(ind_varname), data = subset(df_results, distr=='F')))
  summary(aov(sample_sd ~ get(ind_varname), data = subset(df_results, distr=='F')))
  
  # Statistics of groups
  
  # Compute adjacent pairwsie comparisons, startign with pairwise and testing only adjacent ones
  pw_sample_mean <- pairwise.t.test(subset(df_results, distr=='F')$sample_mean, 
                                    subset(df_results, distr=='F')[[ind_varname]], p.adj = "bonf")$p.value <0.05
  #upperTriangle(pw_sample_mean) = lowerTriangle(pw_sample_mean, byrow=TRUE)
  
  pw_sample_sd <- pairwise.t.test(subset(df_results, distr=='F')$sample_sd, 
                                  subset(df_results, distr=='F')[[ind_varname]], p.adj = "bonf")$p.value <0.05
  #upperTriangle(pw_sample_sd) = lowerTriangle(pw_sample_sd, byrow=TRUE)
  
  
  # Initiialize output struct to store labels of significant for plots
  df_out <- tibble(prox_mean_sig_str = rep("",length(unique(df_results[[ind_varname]]))),
                   prox_sd_sig_str = rep("",length(unique(df_results[[ind_varname]]))),
                   prepost_mean_sig_str = rep("",length(unique(df_results[[ind_varname]]))),
                   prepost_sd_sig_str = rep("",length(unique(df_results[[ind_varname]]))))
  

  # Calculate adjacent significance
  df_out$prox_mean_sig_str <- adjacent_compare_str(pw_sample_mean,'#')
  df_out$prox_sd_sig_str <- adjacent_compare_str(pw_sample_sd,'#')

  # Perform paired ttest between groups at each of the values for swept indepdent variable
  ind_var_levels <- unique(df_results[[ind_varname]])
  p_values_mean <- rep(0,length(ind_var_levels))
  p_values_sd <- rep(0,length(ind_var_levels))
  
  for (n in seq(1,length(ind_var_levels),1)) {
    # Get matching data from 
    fnorm_data <- subset(df_results, df_results[[ind_varname]]==ind_var_levels[n] & df_results$distr=='F')
    norm_data <- subset(df_results, df_results[[ind_varname]]==ind_var_levels[n] & df_results$distr=='N')
    merged_data <- merge(fnorm_data,norm_data, by=c("index",ind_varname))
    
    # Paired ttest (same distribution draws for transform)
    # P values uncorrected at this point
    p_values_mean[n] = t.test(merged_data$sample_mean.x, merged_data$sample_mean.y)$p.value
    p_values_sd[n] = t.test(merged_data$sample_sd.x, merged_data$sample_sd.y)$p.value
    
  }
  # Add chars of significance to groups that are siginificantly different
  df_out$prepost_mean_sig_str = rep('*', length(ind_var_levels))
  df_out$prepost_mean_sig_str[(p.adjust(p_values_mean, method = "bonferroni") > 0.05)] <- ""
  
  df_out$prepost_sd_sig_str = rep('*',length(ind_var_levels))
  df_out$prepost_sd_sig_str[(p.adjust(p_values_sd, method = "bonferroni") > 0.05)] <- ""
  
  # browser();
  
  return(df_out)
  
}


adjacent_compare_str <- function (pw_mat,sig_char) {
#  Given pairwsie comparison table, look for each study group whether it has siginifance with the 
# group before and after
  
  # browser();
  bw_adj_sig = logical(length = dim(pw_mat)[1])
  
  # bw_adj_sig[1] <- pw_mat[1,1]
  
  # bw_adj_sig[length(bw_adj_sig)] <- pw_mat[dim(pw_mat)[1],dim(pw_mat)[2]]
  
  
  
  for (n in seq(1,dim(pw_mat)[1],1)) {
    bw_adj_sig[n] <- pw_mat[n,n] #pw_mat[n-1,n-1] | pw_mat[n,n]
  }
  
  # bw_adj_sig <- pw_mat[,1]
  bw_adj_sig_str = c(rep((sig_char),dim(pw_mat)[1]),"")
  bw_adj_sig_str[!bw_adj_sig] <- ""
  
  # browser();
  return(bw_adj_sig_str)
}
