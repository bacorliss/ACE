



gglabel_height <- function (var,line_num, pad_mult) {
  line_ypos <- pad_mult*(max(var)-min(var)) + max(var) + line_num * pad_mult*(max(var)-min(var))
}




rnorm_swept_param <- function (pop_mu,pop_sigma,n_samples=1e4, n_obs=100) {
# Given mu and sigma, one of them being a number and the other being a vector of values, generate a series of
# simulations of a null distribution, one for each of the elements in the vectored parameter. 
# The other variable remains fixed.
#  
# Inputs
#  mu: population mean for normal distribution, can be single number or vector
#  sigma: population sd for normal distribution, can be a single number or vector
#  n_samples: the number of samples drawn for each simulation
#  n_obs: the number of observations per sample (N)
#  n_sims: number of simulations, equal to the length of the vectored parameter.
  
  # Simulation results to be calculated
  #  Rows: simulation #
  #  Cols: swept mu value
  
# Determine which parameter is vectored to sweep simulations over
param_ind <-1
if (length(pop_mu)<length(pop_sigma)) {param_ind <-2}
param_names <- c("pop_mu","pop_sigma")
  
# Ether sigma or mu are swept
n_sims <- length(get(param_names[param_ind]))

# Simulation results to be calculated
# Rows: simulation #
# Cols: swept mu value
samples_mu <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)
samples_sd <- matrix(rep(0,n_samples*n_sims), ncol = n_sims)


for (n in seq(1,n_sims,1)) {
  # Generate  samples for this simulation, calculate mean and sd for each, store in column
  obs_data <- abs(matrix(rnorm(n_samples*n_obs, mean = pop_mu[min(n,length(pop_mu))], 
                               sd = pop_sigma[min(n,length(pop_sigma))]),
                         nrow = n_trials, byrow = TRUE))
  # Mean of samples for each trial
  samples_mu[,n] <- t(apply(obs_data,1,mean))
  # Standard deviation of samples for each trial
  samples_sd[,n] <- t(apply(obs_data,1,sd))
}

# Collapse matrix of mu's and sd's into a data frame
tbl_mean <-as.table(samples_mu,keep.rownames=FALSE,row.names=NULL,responseName=param_names[param_ind])
colnames(tbl_mean) <- factor(get(param_names[param_ind]))
df_mean <- as.data.frame(tbl_mean) %>% 
  rename(index=Var1,!!eval(param_names[param_ind]) := Var2, "sample_mean" = Freq)

tbl_sd <-as.table(samples_sd, keep.rownames = FALSE,row.names=NULL,responseName=y)
colnames(tbl_sd) <- factor(get(param_names[param_ind]))
df_sd <- as.data.frame(tbl_sd) %>% 
  rename(index=Var1,!!param_names[param_ind] := Var2, sample_sd = Freq)


# Merge mean and sd datatables 
df_results <-merge(df_mean, df_sd, by=c("index",param_names[param_ind]))


}