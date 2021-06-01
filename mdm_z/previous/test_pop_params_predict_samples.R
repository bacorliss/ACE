





sum((df$mu_a + qnorm(0.95) * df$sigma_a/sqrt(df$n_a)) < rowMeans(x_a) )/df$n_samples

sum((df$mu_b + qnorm(0.95) * df$sigma_b/sqrt(df$n_b)) < rowMeans(x_b) )/df$n_samples


sum(df$mu_dm + qnorm(0.95) * df$sigma_dm  < rowMeans(x_b)-rowMeans(x_a) )/df$n_samples
    
    
sum(df_init$rmu_eoc_ucl  < (rowMeans(x_b) / rowMeans(x_a)) )/df$n_samples    
    
    
    
sum(df_init$rmu_dm_ucl  < (rowMeans(x_b) - rowMeans(x_a))/rowMeans(x_a) ) / df$n_samples 
    





# Positive negative conundrum
    
