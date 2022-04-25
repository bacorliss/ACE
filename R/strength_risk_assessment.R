

#' Strength Contests
#' To probe how effective candidate statistics are at quantifying null or effect strength 
#' between study group means. A series of population parameter sets for two 
#' experiments (Exp 1 and 2) with two groups (control group A & experiment 
#' group B) are generated, each a row in the data frame of generated experiment 
#' data. Repeated samples are drawn from the population parameter sets and candidate
#' statistics are quantified and used to determine whether exp 1 or exp 2 has 
#' higher strength
#' All 6 strength parameters (mu_DM, sigma_D, rmu_dm, rsigma_dm, df_D) can be varied as 
#' independent variables to determine how effective candidate statistics are in
#' quantifying strength



# Load required packages
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# Load package manager
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(broom)
p_load(tidyr)
p_load(cowplot)
p_load(dplyr)
# Parallel processing
p_load(boot)
p_load(foreach)
p_load(doParallel)
p_load(stringr)
p_load(confintr)
# User defined libraries
source("R/parallel_utils.R")
source("R/row_stats_toolbox.R")
source("R/aces.R")

# Parse all functions in file for parallel processing using user functions
row_stats_toolbox_fun <- parse_functions_source("R/row_stats_toolbox.R")
aces_functions <- parse_functions_source("R/aces.R")


# Default distribution for population parameters for Exp 1 {a,b}, Exp 2 {a,b}
generate_population_configs <- function(n_samples, n_sims, rand.seed,
                                    # Control group pop. parameters
                                    mus_1a, sigmas_1a, 
                                    mus_2a = NA, sigmas_2a = NA,
                                    # Experiment group pop. parameters
                                    mus_1b = NA, sigmas_1b = NA, 
                                    mus_2b = NA, sigmas_2b = NA,
                                    # Difference distribution pop. parameters
                                    mus_1ao = NA, sigmas_1ao = NA, 
                                    mus_2ao = NA, sigmas_2ao = NA,
                                    # Relative pop. params
                                    rmus_1d = NA, rsigmas_1d = NA,
                                    rmus_2d = NA, rsigmas_2d = NA,
                                    
                                    n_1a = NA, n_1b = NA, n_2a = NA, n_2b = NA,
                                    alpha_1 = 0.05, alpha_2 = 0.05,
                            
                                    toggle_sign_rmu_d_hold_rsigma = FALSE,
                                    toggle_sign_mean_ab = FALSE,
                                    switch_group_ab = FALSE,
                                    switch_mu_ab_12 = FALSE,
                                    switch_mu_d_12 = FALSE,
                                    switch_rmu_d_12_hold_rsigma = FALSE,
                                    switch_sigma_ab_12 = FALSE,
                                    switch_rsigma_ab_12_hold_sigma_a = FALSE,
                                    switch_alpha_12 = FALSE,
                                    switch_n_12 = FALSE,
                                    
                                    fig_name = "test.tiff",
                                    fig_path = "Figure/",
                                    strength_type = "hnst",
                                    gt_colnames, is_plotted = TRUE,
                                    tol = 1e-15) {
  #' @description Generate simulated experiment data for two experiments with 
  #' distributions for mu and sigma specified as vectors (across samples)
  #' 
  #' @param n_samples number of samples drawn per simulation
  #' @param n_sims number of simulations (rows in dataframe)
  #' @param rand.seed seed for random number generation
  #' @param mus_1a means for group a experiment 1
  #' @param sigmas_1a stds for group a experiment 1 
  #' @param mus_2a means for group a experiment 2
  #' @param sigmas_2a stds for group a experiment 2
  #' @param mus_1b means for group b experiment 1
  #' @param sigmas_1b stds for group b experiment 1 
  #' @param mus_2b means for group b experiment 2
  #' @param sigmas_2b stds for group b experiment 2
  #' @param n_1a sample size group a experiment 1
  #' @param n_1b sample size group b experiment 1
  #' @param n_2a sample size group a experiment 2
  #' @param n_2b sample size group b experiment 2
  #' @param mus_1ao mu offset between group a and b experiment 1
  #' @param sigmas_1ao sigma offset between group a and b experiment 1
  #' @param mus_2ao mu offset between group a and b experiment 2
  #' @param sigmas_2ao sigma offset between group a and b experiment 2
  #' @param switch_group_ab flag to randomly switch group a and b assigments
  #' @param toggle_sign_rmu_d_hold_rsigma flag to randomly switch direction of change from a to b
  #' @param switch_exp_12 flag to randomly switch experiment 1 and 2 assignments
  #' @param fig_name base string name of output figures saved to disk
  #' @param fig_path path to output figures saved to disk
  #' @param gt_colnames list of columns in df that serves as independent variables
  #'  and serve as groundtruth for determining whether exp 1 or 2 has higher 
  #'  strength for each pop. param set.
  #' @param is_plotted flag to export figures about pop. param sets to disk
  #' @return df dataframe with generated pop. param sets

  save(list = ls(all.names = TRUE), file = "temp/generate_population_configs.RData",
       envir = environment())
  # load(file = "temp/generate_population_configs.RData")
  
  # Expand any singleton pop param arguments replicate to number of simulations
  input_args <- formalArgs(generate_population_configs)
  pargs <- grep("(^mus)|(^sigmas)|(^alpha)|(^n_\\d)", input_args, value=TRUE)
  # For any pop param not equal in length to n_sims, expand
  for (n in seq_along(pargs)) {
    if (length(get(pargs[n]))==1) assign(pargs[n], rep(get(pargs[n]),n_sims))
  }

  # Record some parameter values for simulations
  set.seed(rand.seed)
  
  # Check that group 1A specified fully
  if (any(is.na(mus_1a) || is.na(sigmas_1a) || is.na(n_1a))) {
    stop("Group 1A not specified, must include: mus_1a, sigma_1a, n_1a")}

  # # Test if input configurations are only for experiment 1 (single experiment setup)
  # is_single_exp = (is.na(mus_2a) && is.na(sigmas_2a) && is.na(n_2a) )
  
  # If have AO values, fill in B and D
  if ( all(!is.na(mus_1ao) && !is.na(mus_2ao) && !is.na(sigmas_1ao) && !is.na(sigmas_2ao))) 
    {
    # Specified:  mus_xa, sigma_xa, mus_xao, sigma_xao
    # Calculated: mus_xb, sigma_xb, mus_xd, sigma_xd 
    df_init <- 
      pop_configs_from_aoffset( n_samples = n_samples, n_sims= n_sims,
                               mus_1a = mus_1a, sigmas_1a = sigmas_1a,  
                               mus_2a = mus_2a, sigmas_2a = sigmas_2a,
                               mus_1ao = mus_1ao, sigmas_1ao = sigmas_1ao, 
                               mus_2ao = mus_2ao, sigmas_2ao = sigmas_2ao,
                               n_1a = n_1a, n_2a = n_2a, n_1b = n_1b, n_2b = n_2b, 
                               alpha_1 = alpha_1, alpha_2 = alpha_2) 
    
    # If relative offset groups specified, calculate params for B and D
  } else if ( all(!is.na(rmus_1d) && !is.na(rmus_2d) && !is.na(rsigmas_1d) && !is.na(rsigmas_2d)))
   {
    # Specified:  mus_xa, sigma_xa, rmus_xd, rsigma_xd
    # Calculated: mus_xb, sigma_xb, mus_xd, sigma_xd 
    
    df_init <- 
      pop_configs_from_rmu_sigma( n_samples = n_samples, n_sims= n_sims,
                               mus_1a = mus_1a, sigmas_1a = sigmas_1a,  
                               mus_2a = mus_2a, sigmas_2a = sigmas_2a,
                               rmus_1d = rmus_1d, rsigmas_1d = rsigmas_1d, 
                               rmus_2d = rmus_2d, rsigmas_2d = rsigmas_2d,
                               n_1a = n_1a, n_2a = n_2a, n_1b = n_1b, n_2b = n_2b, 
                               alpha_1 = alpha_1, alpha_2 = alpha_2) 
    
    
  } else {
    # Specified: mus_xa, sigma_xa, mus_xb, sigma_x_b
    # Calculated: mus_xd, sigma_xd
    # Single experiments (only experiment 1) are processed here
    df_init   <- 
      pop_configs_from_ab( n_samples = n_samples, n_sims = n_sims, 
                          mus_1a = mus_1a, sigmas_1a = sigmas_1a,  mus_2a = mus_2a, sigmas_2a = sigmas_2a,
                          mus_1b = mus_1b, sigmas_1b = sigmas_1b, mus_2b = mus_2b, sigmas_2b = sigmas_2b,
                          n_1a = n_1a, n_2a = n_2a, n_1b = n_1b, n_2b = n_2b, 
                          alpha_1 = alpha_1, alpha_2 = alpha_2) 
  }
  
  # Switch params between groups/experiments if specified with flags
  df <- pop_configs_switches(df = df_init, toggle_sign_rmu_d_hold_rsigma = toggle_sign_rmu_d_hold_rsigma, 
                            toggle_sign_mean_ab = toggle_sign_mean_ab, 
                            switch_group_ab = switch_group_ab,
                            switch_mu_ab_12 = switch_mu_ab_12, 
                            switch_mu_d_12 = switch_mu_d_12,
                            switch_rmu_d_12_hold_rsigma = switch_rmu_d_12_hold_rsigma,
                            switch_sigma_ab_12 = switch_sigma_ab_12, 
                            switch_rsigma_ab_12_hold_sigma_a = switch_rsigma_ab_12_hold_sigma_a,
                            switch_alpha_12 = switch_alpha_12,
                            switch_n_12 = switch_n_12)
  
  # Dataframes that store direction of inequality for each parameter
  # hat: Exp 1 higher null strength than exp 2
  df_hnst = data.frame(is_mudm_1hnst2 = 0)
  # hdt: Exp 2 higher effect strength than exp 2
  df_hest = data.frame(is_mudm_1hest2 = 0)
  
  # Mean of the difference
  
  # Mean and std of difference in means (taken from mean of D since we had the option
  # to invert the sign for D earlier in code)
  df$mu_1dm <- df$mu_1d
  df$mu_2dm <- df$mu_2d
  df$is_mudm_1hnst2 <-  abs(df$mu_1dm) < abs(df$mu_2dm)
  df_hnst$is_mudm_1hnst2 <- "lt"
  df$is_mudm_1hest2 <-  !df$is_mudm_1hnst2
  df_hest$is_mudm_1hest2 <- "gt"
  
  # df$mu_1d <- df$mu_1b - df$mu_1a
  # df$mu_2d <- df$mu_2b - df$mu_2a
  # Is: Exp2 mu[d] > Exp1 mu[d]
  # df$is_mudm_1hnst2 <-  abs(df$mu_1d) < abs(df$mu_2d)
  # df_hnst$is_mudm_1hnst2 <- "lt"
  # df$is_mudm_1hest2 <- !df$is_mudm_1hnst2
  # df_hest$is_mudm_1hest2 <- "gt"
  
  # STD of the difference
  df$sigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2)
  df$sigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) 
  df$is_sigmad_1hnst2 <-  df$sigma_1d < df$sigma_2d
  df_hnst$is_sigmad_1hnst2 <- "lt"
  df$is_sigmad_1hest2 <- df$is_sigmad_1hnst2
  df_hest$is_sigmad_1hest2 <- "lt"
  
  # Degrees of freedom of the difference
  df$df_1d <- df$n_1a + df$n_1b - 2
  df$df_2d <- df$n_2a + df$n_2b - 2
  df$is_dfd_1hnst2 <- df$df_1d > df$df_2d
  df_hnst$is_dfd_1hnst2 <- "gt"
  df$is_dfd_1hest2 <- df$is_dfd_1hnst2
  df_hest$is_dfd_1hest2 <- "gt"
  
  # Degrees of freedom of the difference in means
  df$is_dfdm_1hnst2 <- df$df_1d > df$df_2d
  df_hnst$is_dfdm_1hnst2 <- "gt"
  df$is_dfdm_1hest2 <- df$is_dfdm_1hnst2
  df_hest$is_dfdm_1hest2 <- "gt"
  
  # Pooled standard deviation
  df$sigma_1pool <- sqrt( ( (df$n_1a-1)*df$sigma_1a^2 + (df$n_1b-1)*df$sigma_1b^2) /
                            (df$n_1a-1 + df$n_1b -1 )  )
  df$sigma_2pool <- sqrt( ( (df$n_2a-1)*df$sigma_2a^2 + (df$n_2b-1)*df$sigma_2b^2) /
                            (df$n_2a-1 + df$n_2b-1 ))
  df$is_sigmapool_1hnst2 <- df$sigma_1pool < df$sigma_2pool
  df_hnst$is_sigmapool_1hnst2 <- "lt"
  df$is_sigmapool_1hest2 <-  df$is_sigmapool_1hnst2
  df_hest$is_sigmapool_1hest2 <- "lt"
  
  
  # STD of the difference in means
  df$sigma_1dm <- sqrt(df$sigma_1a^2/df$n_1a + df$sigma_1b^2/df$n_1b)
  df$sigma_2dm <- sqrt(df$sigma_2a^2/df$n_2a + df$sigma_2b^2/df$n_2b)
  df$is_sigmadm_1hnst2 <-  df$sigma_1dm < df$sigma_2dm
  df_hnst$is_sigmadm_1hnst2 <- "lt"
  df$is_sigmadm_1hest2 <- df$is_sigmadm_1hnst2
  df_hest$is_sigmadm_1hest2 <- "lt"
  
  # ALPHA: significance level
  # Higher significance level, higher null strength
  df$is_alpha_1hnst2     <- df$alpha_1 > df$alpha_2
  df_hnst$is_alpha_1hnst2 <- "gt"
  # Higher significance level, higher effect strength
  df$is_alpha_1hest2     <- df$alpha_1 > df$alpha_2
  df_hest$is_alpha_1hest2 <- "gt" 
  
  # Critical t value of difference in means
  df$tstat_1dm <- df$mu_1dm / df$sigma_1dm
  df$tstat_2dm <- df$mu_2dm / df$sigma_2dm
  df$is_tstatdm_1hnst2 <- df$tstat_1dm < df$tstat_2dm
  df_hnst$is_tstatdm_1hnst2 <- "lt"
  df$is_tstatdm_1hest2 <- df$tstat_1dm > df$tstat_2dm
  df_hest$is_tstatdm_1hest2 <- "gt"
  # Critical t value for population parameter set
  df$tcrit_1dm <- qt(1-df$alpha_1, df$n_1a + df$n_1b - 2, lower.tail = TRUE)
  df$tcrit_2dm <- qt(1-df$alpha_2, df$n_2a + df$n_2b - 2, lower.tail = TRUE)
  # T ratio is t statistic dividied by critical t
  df$tratio_1dm <- df$tstat_1dm / df$tcrit_1dm
  df$tratio_2dm <- df$tstat_2dm / df$tcrit_2dm
 

  # Calculate ratio of sigma_md/mu_md to determine how close DM is close to zero,
  # determines whether results are in null region of critical region of t-test
  df$mu_ov_sigma_1dm <- df$mu_1dm / df$sigma_1dm
  df$mu_ov_sigma_2dm <- df$mu_2dm / df$sigma_2dm
  
  # Statistics of difference of means distribution 
  df$rmu_1dm <- df$mu_1dm / df$mu_1a
  df$rmu_2dm <- df$mu_2dm / df$mu_2a
  df$is_rmudm_1hnst2 <-  abs(df$rmu_1dm) - abs(df$rmu_2dm) + tol < 0
  df_hnst$is_rmudm_1hnst2 <- "lt"
  df$is_rmudm_1hest2 <-  abs(df$rmu_1dm) - abs(df$rmu_2dm) - tol > 0
  df_hest$is_rmudm_1hest2 <- "gt"
  
  # Relative sigma of difference
  # browser()
  df$rsigma_1d <- df$sigma_1d / abs(df$mu_1a)
  df$rsigma_2d <- df$sigma_2d / abs(df$mu_2a)
  # df$rsigma_1d <- df$sigma_1d / abs(df$mu_1a + df$mu_1d/2)
  # df$rsigma_2d <- df$sigma_2d / abs(df$mu_2a + df$mu_2d/2)
  df$is_rsigmad_1hnst2 <-  df$rsigma_1d - df$rsigma_2d + tol < 0
  df_hnst$is_rsigmad_1hnst2 <- "lt"
  df$is_rsigmad_1hest2 <-  df$rsigma_1d - df$rsigma_2d + tol < 0
  df_hest$is_rsigmad_1hest2 <- "lt"
  
  # Relative Pooled Sigma
  df$rsigma_1pool <- df$sigma_1pool / abs(df$mu_1a + df$mu_1d/2)
  df$rsigma_2pool <- df$sigma_2pool / abs(df$mu_2a + df$mu_2d/2)
  df$is_rsigmapool_1hnst2 <-  df$rsigma_1pool - df$rsigma_2pool + tol < 0
  df_hnst$is_rsigmapool_1hnst2 <- "lt"
  df$is_rsigmapool_1hest2 <-  df$rsigma_1pool - df$rsigma_2pool + tol < 0
  df_hest$is_rsigmapool_1hest2 <- "lt"
  
  # sigma of the difference of means distribution
  df$rsigma_1dm <- df$sigma_1dm / abs(df$mu_1a)
  df$rsigma_2dm <- df$sigma_2dm / abs(df$mu_2a)
  df$is_rsigmadm_1hnst2 <- df$rsigma_1dm - df$rsigma_2dm + tol < 0
  df_hnst$is_rsigmadm_1hnst2 <- "lt"
  df$is_rsigmadm_1hest2 <- df$rsigma_1dm - df$rsigma_2dm + tol < 0
  df_hest$is_rsigmadm_1hest2 <- "lt"
  
  # Population parameter differences
  df$mean_mud_2m1 <- df$mu_2dm - df$mu_1dm
  df$mean_rmud_2m1 <- (df$mu_2dm/df$mu_2a) - (df$mu_1dm/df$mu_1a)
  df$mean_sigmadm_2m1 <- df$sigma_2dm - df$sigma_1dm
  df$mean_rsigmadm_2m1 <- df$sigma_2dm/df$mu_2a - df$sigma_1dm/df$mu_1a
  
  # Caclulate Power
  # df$power_1 <- pwr.t2n.test(n1 = df$n_1a, n2 = df$n_1a, d = df$mu_1d/df$sigma_1dm, 
  #                            sig.level = df$alpha_1, power = NULL,
  #                            alternative = "two.sided")$power
  # df$power_2 <- pwr.t2n.test(n1 = df$n_2a, n2 = df$n_2a, d = df$mu_2d/df$sigma_2dm, 
  #                            sig.level = df$alpha_2, power = NULL,
  #                            alternative = "two.sided")$power
  
  attr(df,"df_hnst") <- df_hnst
  attr(df,"df_hest") <- df_hest
  
  # browser();
  # Plot generated population parameters
  if (is_plotted){
    plot_population_configs(df, fig_name = fig_name, fig_path = fig_path, 
                           gt_colnames = gt_colnames, strength_type = strength_type)
  } else {
    # browser();
    
    # Plot values of of mu, rmu, sigma, rsigma of d and b over simulations, and df
    df_runs = tibble(Series = rep(seq(1,dim(df)[1],1),6),
                     param = c(rep("mu[DM]",dim(df)[1]), rep("sigma[D]",dim(df)[1]),
                               rep("r*mu[DM]",dim(df)[1]), rep("r*sigma[D]",dim(df)[1]),
                               rep("df[D]",dim(df)[1]),rep("alpha[DM]",dim(df)[1])), 
                     Value = c(df$mu_1dm, df$sigma_1d, df$rmu_1dm, df$rsigma_1d, df$df_1d,df$alpha_1)
    )
    df_runs$param <- factor(df_runs$param, levels = c("mu[DM]", "sigma[D]", "r*mu[DM]","r*sigma[D]","df[D]","alpha[DM]"))
    # Plot facet of each strength measure over series
    gg <- ggplot(data = df_runs, aes(x = Series, y = Value)) +
      geom_line() +
      facet_wrap(vars(param), nrow=3,ncol=2,scales="free_y",labeller=label_parsed) +
      theme_classic(base_size = 8) + ylab("Value") +
      theme(strip.text.x = element_text( margin = margin( b = 0.5, t = 0.5) ),
            panel.spacing = unit(.25, "lines"),)
            # panel.grid.major = element_blank(),
            # panel.grid.minor = element_blank(),
            # strip.background = element_blank(),
            # panel.border = element_rect(fill = NA,colour = "black")) 
    print(gg)
    save_plot(paste(fig_path,  '/', str_replace(fig_name,"\\.[a-z]*$","_params.tiff"), sep = ""), gg, ncol = 1, nrow = 1, 
              base_height = 1.8, base_asp = 4, base_width = 3, dpi = 600)
  }
  return(df)
}

pop_configs_from_aoffset <-
  function( n_samples, n_sims, mus_1a, sigmas_1a, mus_2a, sigmas_2a,
            mus_1ao, sigmas_1ao, mus_2ao, sigmas_2ao,
            n_1a, n_1b, n_2a, n_2b, alpha_1, alpha_2) {
    #' @description Given pop params for group A and offset values (AO), calculates pop params for
    #' group b. This is useful for generating population params because the 
    #' distribution of the offset between group params can be easily 
    #' controlled by explicitly specifying them.
    #' 
    #' mus_ao = mus_b - mus_a
    #' sigmas_ao = sigmas_b - sigmas_a
    #' 
    #' @param n_samples number of samples drawn per simulation
    #' @param n_sims number of simulations (rows in dataframe)
    #' @param mus_1a means for group a experiment 1
    #' @param sigmas_1a stds for group a experiment 1 
    #' @param mus_2a means for group a experiment 2
    #' @param sigmas_2a stds for group a experiment 2
    #' @param mus_1ao mu offset between group a and b experiment 1
    #' @param sigmas_1ao sigma offset between group a and b experiment 1
    #' @param mus_2ao mu offset between group a and b experiment 2
    #' @param sigmas_2ao sigma offset between group a and b experiment 2
    #' @param n_1a sample size group a experiment 1
    #' @param n_1b sample size group b experiment 1
    #' @param n_2a sample size group a experiment 2
    #' @param n_2b sample size group b experiment 2
    #' 
    #' @return df dataframe of population params, where each row is a set. Repeated 
    #' samples are taken with these params and statistics quantified
    
    # Calculate D based on A and offset from A parameters
    mus_1d = mus_1ao;  sigmas_1d = sigmas_1a + sigmas_1ao
    mus_2d = mus_2ao;  sigmas_2d = sigmas_2a + sigmas_2ao
    # calculate B from A and D
    mus_1b = mus_1a + mus_1d;  sigmas_1b = sqrt(sigmas_1d^2 - sigmas_1a^2)
    mus_2b = mus_2a + mus_2d;  sigmas_2b = sqrt(sigmas_2d^2 - sigmas_2a^2)
    # Initialize df with mu_1a, sigma_1a, mu_1b, sigma_1b, mu_1d,
    df = tibble( n_obs = n_obs, n_samples = n_samples,
                 mu_1a = mus_1a, mu_1b = mus_1b, mu_1d = mus_1d,  n_1a = n_1a, n_1b = n_1b,
                 mu_2a = mus_2a, mu_2b = mus_2b, mu_2d = mus_2d,  n_2a = n_2a, n_2b = n_2b, 
                 sigma_1a = sigmas_1a, sigma_1b = sigmas_1b, sigma_1d = sigmas_1d,
                 sigma_2a = sigmas_2a, sigma_2b = sigmas_2b, sigma_2d = sigmas_2d,
                 alpha_1 = alpha_1, alpha_2 = alpha_2 
    )
    return(df)
  }



pop_configs_from_rmu_sigma <-
  function( n_samples, n_sims, mus_1a, sigmas_1a, mus_2a, sigmas_2a,
            rmus_1d, rsigmas_1d, rmus_2d, rsigmas_2d,
            n_1a, n_1b, n_2a, n_2b, alpha_1, alpha_2) {
    
    save(list = ls(all.names = TRUE), file = "temp/pop_configs_from_rmu_sigma.RData",envir = environment())
    # load(file = "temp/pop_configs_from_rmu_sigma.RData")
    
    
    # Calculate mu and sigma d and b based on relative mus and sigmas
    mus_1d = rmus_1d * mus_1a;  mus_2d = rmus_2d * mus_2a; 
    # Calculate mu for group b
    mus_1b =  mus_1d + mus_1a;  mus_2b =  mus_2d + mus_2a;
      
    sigmas_1d = rsigmas_1d * mus_1a 
    sigmas_2d = rsigmas_2d * mus_2a 
    
    sigmas_1b = sqrt(sigmas_1d^2 - sigmas_1a^2)
    sigmas_2b = sqrt(sigmas_2d^2 - sigmas_2a^2)
    
    # Initialize df with mu_1a, sigma_1a, mu_1b, sigma_1b, mu_1d,
    df = tibble( n_obs = n_obs, n_samples = n_samples,
                 mu_1a = mus_1a, mu_1b = mus_1b, mu_1d = mus_1d,  n_1a = n_1a, n_1b = n_1b,
                 mu_2a = mus_2a, mu_2b = mus_2b, mu_2d = mus_2d,  n_2a = n_2a, n_2b = n_2b, 
                 sigma_1a = sigmas_1a, sigma_1b = sigmas_1b, sigma_1d = sigmas_1d,
                 sigma_2a = sigmas_2a, sigma_2b = sigmas_2b, sigma_2d = sigmas_2d,
                 alpha_1 = alpha_1, alpha_2 = alpha_2
    )
  }





pop_configs_from_ab <- function( n_samples, n_sims, 
                                mus_1a, sigmas_1a, 
                                mus_2a, sigmas_2a,
                                mus_1b, sigmas_1b, 
                                mus_2b, sigmas_2b,
                                n_1a, n_1b, n_2a, n_2b, alpha_1, alpha_2) {
  #' @description Given pop params for group A and B, calculates pop params of D, 
  #' the difference between A and B distribution.
  #' 
  #' mus_d = mus_b - mus_a
  #' sigmas_d = sqrt(sigmas_a^2 + sigmas^2)
  #' 
  #' @param n_samples number of samples drawn per simulation
  #' @param n_sims number of simulations (rows in dataframe)
  #' @param mus_1a means for group a experiment 1
  #' @param sigmas_1a stds for group a experiment 1 
  #' @param mus_2a means for group a experiment 2
  #' @param sigmas_2a stds for group a experiment 2
  #' @param mus_1b means for group b experiment 1
  #' @param sigmas_1b stds for group b experiment 1 
  #' @param mus_2b means for group b experiment 2
  #' @param sigmas_2b stds for group b experiment 2
  #' @param n_1a sample size group a experiment 1
  #' @param n_1b sample size group b experiment 1
  #' @param n_2a sample size group a experiment 2
  #' @param n_2b sample size group b experiment 2
  #' 
  #' @return df dataframe of population params, where each row is a set.
  
  # Calculate D based on A and B parameters
  mus_1d = mus_1b - mus_1a;  sigmas_1d = sqrt(sigmas_1a^2 + sigmas_1b^2)
  mus_2d = mus_2b - mus_2a;  sigmas_2d = sqrt(sigmas_2a^2 + sigmas_2b^2)
  # Initialize df with mu_1a, sigma_1a, mu_1b, sigma_1b, mu_1d,
  df = tibble( n_obs = n_obs, n_samples = n_samples,
               mu_1a = mus_1a, mu_1b = mus_1b, mu_1d = mus_1d,  n_1a = n_1a, n_1b = n_1b,
               mu_2a = mus_2a, mu_2b = mus_2b, mu_2d = mus_2d,  n_2a = n_2a, n_2b = n_2b, 
               sigma_1a = sigmas_1a, sigma_1b = sigmas_1b, sigma_1d = sigmas_1d,
               sigma_2a = sigmas_2a, sigma_2b = sigmas_2b, sigma_2d = sigmas_2d,
               alpha_1 = alpha_1, alpha_2 = alpha_2
  )
  return(df)
}


pop_configs_switches <- function(df_init, toggle_sign_rmu_d_hold_rsigma, toggle_sign_mean_ab, 
                                switch_group_ab, switch_mu_ab_12, switch_rmu_d_12_hold_rsigma,
                                switch_mu_d_12, switch_sigma_ab_12, switch_rsigma_ab_12_hold_sigma_a,
                                switch_alpha_12, switch_n_12) {
  #' @description Apply various switches/alterations to group assignment of 
  #' pop. param sets at random. These alterations allow the strength measures 
  #' to be probed in an independent fashion by averaging out the effects of all 
  #' other strength measures
  #' 
  #' @param df_init dataframe of population params
  #' @param toggle_sign_rmu_d_hold_rsigma flag to randomly switch the sign of the difference 
  #' in means between A and B. So if A is larger than a by some amount, B now 
  #' becomes smaller than A by the same amount
  #' @param toggle_sign_mean_ab flag to randomly switch the sign of the means 
  #' for A and B
  #' @param switch_group_ab flag to randomly switch group assignments between 
  #' A and B 
  #' @param switch_exp_12 means for group a experiment 2
  #' 
  #' @return df dataframe of pop param sets with switches applied 
  
  df <- df_init
  
  if (switch_rsigma_ab_12_hold_sigma_a) {
    # browser();
    # Switch sigma_a, sigma_b, sigma_d between experiments by switching rsigma_a, rsigma_b
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    # Switch sigma_a and sigma_b
    prev_sigma_1a <- df$sigma_1a;   
    prev_sigma_2a <- df$sigma_2a;   
    
    prev_rsigma_1d <- df$sigma_1d/df$mu_1a;
    prev_rsigma_2d <- df$sigma_2d/df$mu_2a;
    
    # Reassign sigma_d
    new_sigma_1d <- prev_rsigma_1d * df$mu_2a
    new_sigma_2d <- prev_rsigma_2d * df$mu_1a
    
    
    # Recalculate sigma_b from new rsigma
    df$sigma_1b[switch_boolean]    <- sqrt(new_sigma_2d[switch_boolean]^2 -
                                             prev_sigma_1a[switch_boolean]^2)
    df$sigma_2b[switch_boolean]    <- sqrt(new_sigma_1d[switch_boolean]^2 -
                                             prev_sigma_2a[switch_boolean]^2)
    
    # Recalculate sigma_d from switched rows
    df$sigma_1d[switch_boolean] = new_sigma_2d[switch_boolean]
    df$sigma_2d[switch_boolean] = new_sigma_1d[switch_boolean]
  }
  

  # Randomly switch sign of both group a and b for exp 1 and 2 separately, recalculate d
  if (toggle_sign_mean_ab) {
  # Randomly switch sign of both mu_a, mu_b, and mu_d
    switch_boolean <- sample( c(TRUE,FALSE), n_sims, TRUE)
    prev_mu_1a     <- df$mu_1a; prev_mu_1b     <- df$mu_1b;
    prev_mu_2a     <- df$mu_2a; prev_mu_2b     <- df$mu_2b
    df$mu_1a[switch_boolean] <- -prev_mu_1a[switch_boolean]
    df$mu_2a[switch_boolean] <- -prev_mu_2a[switch_boolean]
    df$mu_1b[switch_boolean] <- -prev_mu_1b[switch_boolean]
    df$mu_2b[switch_boolean] <- -prev_mu_2b[switch_boolean]
    
    # Switch mu_d
    prev_mu_1d <- df$mu_1d;  prev_mu_2d <- df$mu_2d;  
    df$mu_1d[switch_boolean] = -prev_mu_1d[switch_boolean]; 
    df$mu_2d[switch_boolean] = -prev_mu_2d[switch_boolean];
  }
  
  if (switch_group_ab) {
    # Randomly switch mu_a <-> mu_b, and sigma_a <-> sigma_b for each experiment
    # mu_d will also switch signs, while sigma_d will remain unchanged
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    
    # Save temp variables for switch
    prev_mu_1a     <- df$mu_1a;    prev_sigma_1a  <- df$sigma_1a
    prev_mu_1b     <- df$mu_1b;    prev_sigma_1b  <- df$sigma_1b
    prev_mu_2a     <- df$mu_2a;    prev_sigma_2a  <- df$sigma_2a
    prev_mu_2b     <- df$mu_2b;    prev_sigma_2b  <- df$sigma_2b
    # Exp 1
    df$mu_1a[switch_boolean]      <- prev_mu_1b[switch_boolean]
    df$sigma_1a[switch_boolean]   <- prev_sigma_1b[switch_boolean]
    df$mu_1b[switch_boolean]      <- prev_mu_1a[switch_boolean]
    df$sigma_1b[switch_boolean]   <- prev_sigma_1a[switch_boolean]
    # Exp 2
    df$mu_2a[switch_boolean]      <- prev_mu_2b[switch_boolean]
    df$sigma_2a[switch_boolean]   <- prev_sigma_2b[switch_boolean]
    df$mu_2b[switch_boolean]      <- prev_mu_2a[switch_boolean]
    df$sigma_2b[switch_boolean]   <- prev_sigma_2a[switch_boolean]
    
    # Switch mu_d by simply negating sign 
    prev_mu_1d <- df$mu_1d;  prev_mu_2d <- df$mu_2d;  
    df$mu_1d[switch_boolean] = -prev_mu_1d[switch_boolean]
    df$mu_2d[switch_boolean] = -prev_mu_2d[switch_boolean]
    
  }
  if (switch_mu_ab_12) {
    # Switch both mu_a and mu_b and mu_d between experiments
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    
    # Switch mu_a and mu_b
    prev_mu_1a <- df$mu_1a;  prev_mu_1b <- df$mu_1b;  
    prev_mu_2a <- df$mu_2a;  prev_mu_2b <- df$mu_2b;
    df$mu_1a[switch_boolean]      <- prev_mu_2a[switch_boolean]
    df$mu_1b[switch_boolean]      <- prev_mu_2b[switch_boolean]
    df$mu_2a[switch_boolean]      <- prev_mu_1a[switch_boolean]
    df$mu_2b[switch_boolean]      <- prev_mu_1b[switch_boolean]

    # Switch mu_d
    prev_mu_1d <- df$mu_1d;  prev_mu_2d <- df$mu_2d;  
    df$mu_1d[switch_boolean] = prev_mu_2d[switch_boolean]
    df$mu_2d[switch_boolean] = prev_mu_1d[switch_boolean]
  } 
  
  if (switch_mu_d_12) {
    # Switch both mu_a and mu_b and mu_d between experiments
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    # Switch mu_d
    prev_mu_1d <- df$mu_1d;  prev_mu_2d <- df$mu_2d;  
    df$mu_1d[switch_boolean] = prev_mu_2d[switch_boolean]
    df$mu_2d[switch_boolean] = prev_mu_1d[switch_boolean]
    #Recalculate group b mean from d
    df$mu_1b[switch_boolean] <- df$mu_1a[switch_boolean] + df$mu_1d[switch_boolean]
    df$mu_2b[switch_boolean] <- df$mu_2a[switch_boolean] + df$mu_2d[switch_boolean]
  } 
  
  if (switch_sigma_ab_12) {
    # Switch sigma_a, sigma_b, sigma_d between experiments
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    
    # Switch sigma_a and sigma_b
    prev_sigma_1a <- df$sigma_1a;   prev_sigma_1b <- df$sigma_1b;
    prev_sigma_2a <- df$sigma_2a;   prev_sigma_2b <- df$sigma_2b;
    df$sigma_1a[switch_boolean]    <- prev_sigma_2a[switch_boolean]
    df$sigma_1b[switch_boolean]    <- prev_sigma_2b[switch_boolean]
    df$sigma_2a[switch_boolean]    <- prev_sigma_1a[switch_boolean]
    df$sigma_2b[switch_boolean]    <- prev_sigma_1b[switch_boolean]
    
    # Switch sigma_d 
    prev_sigma_1d <- df$sigma_1d;  prev_sigma_2d <- df$sigma_2d;  
    df$sigma_1d[switch_boolean] = prev_sigma_2d[switch_boolean]
    df$sigma_2d[switch_boolean] = prev_sigma_1d[switch_boolean]
  }
  if (switch_alpha_12) {
    # browser();
    prev_alpha_1 <- df$alpha_1;    prev_alpha_2   <-  df$alpha_2;
    # Determine which simulations to switch parameters for a and b
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    df$alpha_1[switch_boolean]    <- prev_alpha_2[switch_boolean]
    df$alpha_2[switch_boolean]    <- prev_alpha_1[switch_boolean]
  }
  if (switch_n_12) {
    # browser();
    prev_n_1a <- df$n_1a;   prev_n_1b <- df$n_1b
    prev_n_2a <- df$n_2a;   prev_n_2b <- df$n_2b
    # Determine which simulations to switch parameters for a and b
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    df$n_1a[switch_boolean]    <- prev_n_2a[switch_boolean]
    df$n_1b[switch_boolean]    <- prev_n_2b[switch_boolean]
    df$n_2a[switch_boolean]    <- prev_n_1a[switch_boolean]
    df$n_2b[switch_boolean]    <- prev_n_1b[switch_boolean]
  }
  
  
  if (switch_rmu_d_12_hold_rsigma) {
    # browser();
    # Calculate mu_d based on switching rmu_d
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    orig_rmu_1d <- df$mu_1d/df$mu_1a; orig_rmu_2d <- df$mu_2d/df$mu_2a;  
   
    
    orig_rsigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2) / df$mu_1a
    orig_rsigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) / df$mu_2a
    

    df$mu_1d[switch_boolean] = orig_rmu_2d[switch_boolean] * df$mu_1a[switch_boolean]
    df$mu_2d[switch_boolean] = orig_rmu_1d[switch_boolean] * df$mu_2a[switch_boolean]
    #Recalculate group b mean from d
    df$mu_1b[switch_boolean] <- df$mu_1a[switch_boolean] + df$mu_1d[switch_boolean]
    df$mu_2b[switch_boolean] <- df$mu_2a[switch_boolean] + df$mu_2d[switch_boolean]
    
    # Hold rsigma_d constant by calculating sigma_1b
    new_sigma_1d <- orig_rsigma_1d * df$mu_1a
    new_sigma_2d <- orig_rsigma_2d * df$mu_2a
    new_sigma_1b <- sqrt(new_sigma_1d^2 - df$sigma_1a^2)
    new_sigma_2b <- sqrt(new_sigma_2d^2 - df$sigma_2a^2)
    
    new_rsigma_1d <- sqrt(df$sigma_1a^2 + new_sigma_1b^2)/df$mu_1a
    new_rsigma_2d <- sqrt(df$sigma_2a^2 + new_sigma_2b^2)/df$mu_2a
    
    # Change sigma_[1,2]b to that the rsigma is the same
    df$sigma_1b = new_sigma_1b
    df$sigma_2b = new_sigma_2b
    # new_rsigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2) / (df$mu_1a + 0.5 * df$mu_1d)
    # new_rsigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) / (df$mu_2a + 0.5 * df$mu_2d)
  } 
  

  if (toggle_sign_rmu_d_hold_rsigma) { 
    # Randomly switch sign of D for both experiments, recalculate B
    # Then recalculate sigma_b so that rsigma_d is the same as before
    # browser();
    
    swith_boolean = sample(c(TRUE,FALSE), n_sims, TRUE)
    
    orig_mu_1b <- df$mu_1b;  orig_mu_2b <- df$mu_2b;
    orig_mu_1d <- df$mu_1d;  orig_mu_2d <- df$mu_2d;  
    orig_sigma_1b <- df$sigma_1b; orig_sigma_2b <- df$sigma_2b;
    
    orig_rsigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2) / (df$mu_1a)
    orig_rsigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) / (df$mu_2a)
    
    # Switch sign of mu_d randomly
    df$mu_1d[swith_boolean] =  -1 * orig_mu_1d[swith_boolean]
    df$mu_2d[swith_boolean] =  -1 * orig_mu_2d[swith_boolean]
    # Recalculate b based on switched d
    df$mu_1b[swith_boolean] = df$mu_1a[swith_boolean] + df$mu_1d[swith_boolean]
    df$mu_2b[swith_boolean] = df$mu_2a[swith_boolean] + df$mu_2d[swith_boolean]
    
    # Adjust sigma_[1,2]b so that rsigma_1d and rsigma_2d (switched mus only)
    # are constant 
    #
    # rsigma_d =   sigma_d_old     =     sigma_d_new
    #            ______________         ________________
    #              avg_mu_ab_old         avg_mu_ab_new
    #
    # Note: avg_mean_ab = (mu_a + mu_b) / 2
    #
    #
    #  sigma_d_new =  sigma_d_old * avg_mu_ab_new
    #                ____________________________
    #                       avg_mu_ab_old
    #
    #
    #
    
    # rsigma must stay equal, so sigma_b is adjusted
    new_sigma_1d <- orig_rsigma_1d * df$mu_1a
    new_sigma_2d <- orig_rsigma_2d * df$mu_2a
    new_sigma_1b <- sqrt(new_sigma_1d^2 - df$sigma_1a^2)
    new_sigma_2b <- sqrt(new_sigma_2d^2 - df$sigma_2a^2)
    
    new_rsigma_1d <- sqrt(df$sigma_1a^2 + new_sigma_1b^2)/(df$mu_1a)
    new_rsigma_2d <- sqrt(df$sigma_2a^2 + new_sigma_2b^2)/(df$mu_2a)
    
    # Change sigma_[1,2]b to that the rsigma is the same
    df$sigma_1b = new_sigma_1b
    df$sigma_2b = new_sigma_2b
    
  }
  
  # 
  # browser();
  return(df)
}


plot_population_configs <- function(df_init, gt_colnames,fig_name,fig_path, strength_type){
  #' @description plots characterizing selection of pop. params across simulations.
  #' 1) Plot fraction of pop param sets with exp 1 higher strength than 
  #' experiment 2 according to each candidate statistic, and test for non random
  #'  correlation between strength parameters.
  #' 2) Plot histogram of mu[DM]/sigma[DM] to visualize whether selected pop. 
  #' params fall within null or critical region (threshold ~2.5 stds from mean).
  #'
  #' @param df_init initial dataframe of pop params, each row is a set
  #' @param gt_colnames list of columns in df that serves as independent variables
  #'  and serve as groundtruth for determining whether exp 1 or 2 has higher 
  #'  strength for each pop. param set.
  #' @param fig_name
  #' @param fig_path path to output figures saved to disk
  #' 
  #' @return no return, exports figures to disk

  save(list = ls(all.names = TRUE), file = "temp/plot_population_configs.RData",envir = environment())
  # load(file = "temp/plot_population_configs.RData")
  
  # Output csv of strength of input parameters to each individual input parameter
  param_fields = paste(c("is_mudm_1","is_rmudm_1","is_sigmad_1", "is_rsigmad_1",
                         "is_dfdm_1", "is_alpha_1"), strength_type,"2",sep="")


  # bv_gt_colnames <- sapply(gt_colnames, function(x) any(x==param_fields))
  if (!all(sapply(gt_colnames, function(x) any(x==param_fields)))) {stop("gt_colnames is not correctly worded")}
  
  
  # Calculate indices of colname
  gt_param_inds <- unname(sapply(gt_colnames,function(x){pmatch(x,param_fields)}))
  gt_param_labels <- c("abs(~mu[DM]*phantom(.))",  "abs(~r*mu[DM]*phantom(.))",
                       "~sigma[D]*phantom(.)", "r*sigma[D]", "df[D]","alpha[DM]")

  # Calculate agreement matrix: test if nonrandom agreement between each parameter
  # versus every other.
  n_agreement = matrix(0, ncol = length(param_fields), nrow = length(param_fields))
  pwise_binom_p <- n_agreement
  for (r in seq(1,length(param_fields),1)) {
    for (c in seq(1,length(param_fields),1)) {
      n_agreement[r,c]   <- sum(df_init[[param_fields[r]]] == df_init[[param_fields[c]]])
      # Bonferroni corrected binomial proportion test
      pwise_binom_p[r,c] <-
        prop.test(n_agreement[r,c], dim(df_init)[1], alternative = "two.sided",
                  conf.level = 1-0.05/(length(param_fields)*length(gt_colnames)), 
                  correct = TRUE)$p.value
    }
  }
  # Corrected p values based on number of independent variables selected * 4
  pwise_binom_p_corr <- pmin(4*length(gt_colnames) * pwise_binom_p,rep(1,prod(dim(pwise_binom_p))))
  pwise_binom_p_sig <- ifelse(pwise_binom_p_corr<0.05, "#","")
  
  str_binom_p <- matrix(sapply(pwise_binom_p_corr, function(x) if(x> 0.05) 
  {sprintf("%0.2f",x)} else {sprintf("%0.2e",x)}),nrow = dim(pwise_binom_p)[1])
  colnames(str_binom_p) <- str_replace(str_replace(param_fields, "is_", ""), "_2[a-z][a-z]1","")
  rownames(str_binom_p) <- str_replace(str_replace(param_fields, "is_", ""), "_2[a-z][a-z]1","")
  # Write to csv
  csv_path <- paste(fig_path, "params_",str_replace(fig_name, ".tiff$", ".csv"), sep="")
  cat("Table 1: Binomial test of agreement by group\n", file = csv_path)
  suppressWarnings(write.table(str_binom_p, csv_path, append = FALSE, 
                               col.names = TRUE, sep=","))
  
  # Calculate binomial confidence intervals for each success rate for each parameter
  df_params <- tibble(group = c("mu", "rmu", "sigma", "rsigma", "df","alpha"), 
                      estimate = rep(0,length(param_fields)), lci = rep(0,length(param_fields)),
                      uci = rep(0,length(param_fields)))
  for (n in seq_along(param_fields)) {
    binom <- prop.test(sum(df_init[[param_fields[n]]]), dim(df_init)[1], 
                                  conf.level=1-0.05/(5*length(gt_colnames)), correct = FALSE)
    df_params$estimate[n]  <-binom$estimate
    df_params$lci[n]       <-binom$conf.int[1]
    df_params$uci[n]       <-binom$conf.int[2]
  }
  df_params$group <- factor(df_params$group, levels = c("mu", "rmu", "sigma", "rsigma", "df","alpha"))
  
  # Plot confidence interval of success rate for each parameter and their agreement
  gg <- ggplot(df_params, aes(x = group,  y = estimate)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey",linetype="dashed") +
    geom_linerange(aes(ymin = lci, ymax = uci), size = 0.5) +
    geom_point(size = 1.25, fill = "white", shape = 1) + 
    ylab(bquote(Frac.~Exp[1]~.(toupper(strength_type))~Exp[2]~~~~~~~~~~~~phantom("."))) +
    xlab("Groundtruth") +
    theme_classic(base_size = 8) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(y=c(0,1), clip = "off") +
    scale_x_discrete( labels= parse(text = gt_param_labels)) +
    theme(plot.margin = unit(c(13+6.5*(length(gt_colnames)-1),6,3,3), "pt"))

  for (n in seq_along(gt_colnames)) {
    gg <- gg +
      annotate("text",label = paste(gt_param_labels[gt_param_inds[n]],"~phantom(.)",
                                    sep = ""), x = 0.1, size=2.5,
               y = 1+.15*n, vjust = 0.5, hjust = 1, parse = TRUE) +
      geom_text(y = 1+0.15*n, aes(group = n),label = ifelse(pwise_binom_p_corr[gt_param_inds[n],]<0.05, "#",""),
                size = 2.5, vjust=0.5, hjust=0.5) +
      annotate("segment", x = 0.1, xend = length(param_fields)+0.5, y = 1+0.15*n-.08, 
               yend = 1+0.15*n-.08, colour = "black", size=.2) 
    save_plot(paste(fig_path, '/gt_',fig_name, ".tiff", sep = ""), gg, ncol = 1, nrow = 1, 
              base_height = 1.5, base_asp = 3, base_width = 2, dpi = 600)
    
  }
  print(gg)
  save_plot(paste(fig_path, '/gt_',fig_name, ".tiff", sep = ""), gg, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 2, dpi = 600)


  # Export csv file for agreement between each variable to others
  # Plot histogram of mu[D]/sigma[D] to demonstrate how far from zero D is
  df <-tibble(group = as.factor(c(rep(1,dim(df_init)[1]),rep(2,dim(df_init)[1]))),
              t_ratio = c(df_init$tratio_1dm, df_init$tratio_2dm),
              mu_ov_sigma = c(df_init$mu_1dm/df_init$sigma_1dm,
                              df_init$mu_2dm/df_init$sigma_2dm))
  p <- ggplot(df, aes(x = t_ratio, y = t_ratio, fill = group)) +
    geom_histogram(aes(y=stat(count / sum(count))), position = "identity", 
                   alpha=0.25, bins = 30) +
    geom_vline(xintercept = -1,linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    xlab( expression(t[stat]*phantom(.)/phantom(.)*abs(phantom(.)*t[crit]*phantom(.)))) +
    ylab( "Freq.") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01),expand = c(0, 0))
  ymax <- max(ggplot_build(p)$data[[1]]$ymax)
  # KS Test between both groups as they are plotted to see if hist. are different
  fill_id <- levels(as.factor(ggplot_build(p)$data[[1]]$fill))
  ks_p_val <- ks.test(subset(ggplot_build(p)$data[[1]], fill == fill_id[1])$y, 
                      subset(ggplot_build(p)$data[[1]], fill == fill_id[2])$y, 
                      alternative = "two.sided", exact = FALSE)
  p <- p + expand_limits(y = c(0, 1.1* ymax)) + 
    annotate("label",label=ifelse(ks_p_val$p.value>0.05, 
                                  sprintf('p = %.3f', ks_p_val$p.value),
                                  sprintf('p = %.2e', ks_p_val$p.value))
               , x=mean(range(df$t_ratio)), y=1.07*ymax, 
             size=2, fill = "white",label.size = NA)
  # print(p)
  save_plot(paste(fig_path, '/predicted_t-ratio_',fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 1.2, dpi = 600)
  
}


quantify_population_configs <- function(df_in, overwrite = TRUE,
                                out_path = "temp/", data_file_name,
                                rand.seed = 0, 
                                parallel_sims = TRUE,  stat_exclude_list = NULL,
                                strength_type = "hnst", delta, is_delta_relative = FALSE,
                                use_pseudo_samples = use_pseudo_samples) {
  #' @description Given a data frame of pop. params and input data for each 
  #' experiment group, quantify mean values and comparison error of various 
  #' candidate statistics over repeated samples.
  #' 
  #' @param overwrite flag to overwrite data file storing results
  #' @param out_path path to export data file to disk storing results
  #' @param data_file_name filename of exported datafile storing results
  #' @param rand.seed seed for random number generation
  #' @param parallel_sims flag to parallelize simulation across CPU cores (recommended)
  #' 
  #' @return df dataframe with pop params and comparison error rates of each candidate statistics.
  
  # Initialize output data frame
  n_sims = dim(df_in)[1]
  # df <- df_in
  
  save(list = ls(all.names = TRUE), file = "temp/quantify_population_configs.RData",envir = environment())
  # load(file = "temp/quantify_population_configs.RData")

  
  # Only perform simulations if results not saved to disk
  if (!file.exists(paste(out_path,'/',data_file_name,sep="")) | overwrite) {
    
    if (parallel_sims) { #print("starting parallel processing")
      # Setup parallel back end to use many processors
      cl = makeCluster(detectCores()[1]-1)
      registerDoParallel(cl)
      
      df <- foreach(n = 1:n_sims, .combine = rbind, .packages = 
                      c("BayesFactor","TOSTER", "sgpv"),
                .export = c(row_stats_toolbox_fun, aces_functions,
                            "quantify_population_config")) %dopar% {
        #calling a function
        tempMatrix <- quantify_population_config(df_in[n,], rand.seed = rand.seed,
                                               stat_exclude_list = stat_exclude_list,
                                               strength_type = strength_type, delta = delta,
                                               is_delta_relative = is_delta_relative,
                                               use_pseudo_samples = use_pseudo_samples) 
        tempMatrix
      }
      #stop cluster
      stopCluster(cl)
    }else{
      df_list = list();
      # Process effect sizes serially and bind rows into one dataframe
      for (n in seq(1,n_sims)) {
        print(n)
        df_list[[n]] <- quantify_population_config(df_in[n,], rand.seed = rand.seed, 
                                            stat_exclude_list = stat_exclude_list,
                                            strength_type = strength_type, delta = delta,
                                            is_delta_relative = is_delta_relative,
                                            use_pseudo_samples = use_pseudo_samples) 
      }
      df <- bind_rows(df_list)
    }
    # browser();
    
    # Save dataframe results to a file
    saveRDS(df, file = paste(out_path,'/', data_file_name, sep=""))
  } else {
    # Restore the dataframe results from disk
    df <- readRDS(file = paste(out_path,'/', data_file_name,sep=""))
  }
  
  return(df)
}


quantify_population_config <- function(df, rand.seed = 0, 
                                      parallelize_bf = FALSE, stat_exclude_list = NULL, 
                                      strength_type = "hnst", delta, is_delta_relative = FALSE,
                                      use_pseudo_samples = FALSE) {
  #' @description Quantifies mean and comparison error of various candidate 
  #' statistics across repeated samples specified by input population parameters
  #' 
  #' @param df input dataframe with a single row
  #' @param rand.seed seed for random number generation
  #' @param parallelize_bf  flag to parallize calculation of bayes factor 
  #' (not recommended since each pop. param set is already parallelized 
  #' in parent function call).
  #' 
  #' @return df_out output dataframe with a single row, with pop. params and the 
  #' mean and std of each candidate statistic across samples along with their 
  #' comparison error for determining higher strength.
  
  save(list = ls(all.names = TRUE), file = "temp/quantify_population_config.RData",envir = environment())
  # load(file = "temp/quantify_population_config.RData")
  
  if (dim(df)[1] != 1) stop("Need to input a single row of df")
  set.seed(rand.seed)
  
 
  
  # Pseudo sampling
  # Since we are comparing repeated pairs of samples from two experiments, we 
  # can generate a small number of samples from both experiemnts and then compare
  # them in a pairwise fashion. Pseudo_indices are the index from each series of samples
  if (use_pseudo_samples) {
    tot_samples <- df$n_samples
    # Calculate minimum number of n_actual_samples needed to get total samples
    df$n_actual_samples <- ceiling(sqrt(df$n_samples))
    
    all_pseudo_sample_indices <- expand.grid(exp_1 = seq(1,df$n_actual_samples,1), exp_2 = seq(1,df$n_actual_samples,1))
    if (dim(all_pseudo_sample_indices)[1]<df$n_samples) {stop("not enough pseudo samples generated")}
    row_pseudo_indicies <- sample.int(dim(all_pseudo_sample_indices)[1], tot_samples, replace = FALSE)
    pseudo_indicies <- all_pseudo_sample_indices[row_pseudo_indicies,]

  }else{
    tot_samples <- df$n_samples
    df$n_actual_samples <- df$n_samples
    pseudo_indicies = data.frame(exp_1 = seq(1,df$n_actual_samples,1), exp_2 = seq(1,df$n_actual_samples,1))
  }
  

  
  # Generate samples and quantify candidate statistics for experiment 1
  x_1a = matrix(rnorm(df$n_actual_samples * df$n_1a, mean = df$mu_1a, 
                      sd = df$sigma_1a), nrow = df$n_actual_samples, 
                ncol = df$n_1a)
  x_1b = matrix(rnorm(df$n_actual_samples * df$n_1b, mean = df$mu_1b, 
                      sd = df$sigma_1b), nrow = df$n_actual_samples, 
                ncol = df$n_1b)
  dfs_1 <- 
    quantify_row_stats(x_a = x_1a, x_b = x_1b, parallelize_bf = FALSE, 
                       delta = delta, is_delta_relative = is_delta_relative,
                       stat_exclude_list = stat_exclude_list, 
                       conf.level = 1 - df$alpha_1)
  stat_list <- colnames(dfs_1)
  
  
  # Generate samples and quantify candidate statistics for experiment 1
  if (!is.na(df$mu_2a)) {
  x_2a = matrix(rnorm(df$n_actual_samples * df$n_2a, mean = df$mu_2a, 
                      sd = df$sigma_2a), nrow = df$n_actual_samples, 
                ncol = df$n_2a)
  x_2b = matrix(rnorm(df$n_actual_samples * df$n_2b, mean = df$mu_2b, 
                      sd = df$sigma_2b), nrow = df$n_actual_samples, 
                ncol = df$n_2b)
  dfs_2 <- 
    quantify_row_stats(x_a = x_2a, x_b = x_2b, parallelize_bf = FALSE, 
                       delta = delta, is_delta_relative = is_delta_relative,
                       stat_exclude_list = stat_exclude_list, 
                       conf.level = 1 - df$alpha_2)
  } else {
    dfs_2 <- dfs_1
  }
  
  
  
  # save(list = ls(all.names = TRUE), file = "temp/quantify_population_config.RData",envir = environment())
  # load(file = "temp/quantify_population_config.RData")
  
  dfc <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), paste("exp1_mean_",stat_list[1],sep=''))
  # Means and std of each statistic

  for (i in seq_along(stat_list)) {
    dfc[[paste("exp1_mean_", stat_list[i],sep='')]] <- mean(dfs_1[[stat_list[i]]], na.rm = TRUE)
    dfc[[paste("exp1_sd_", stat_list[i],sep='')]]   <- sd(dfs_1[[stat_list[i]]], na.rm = TRUE)
    
    dfc[[paste("exp2_mean_", stat_list[i],sep='')]] <- mean(dfs_2[[stat_list[i]]], na.rm = TRUE)
    dfc[[paste("exp2_sd_", stat_list[i],sep='')]]   <- sd(dfs_2[[stat_list[i]]], na.rm = TRUE)
    
    dfc[[paste("exp12_nincluded_", stat_list[i],sep='')]] <- sum(!is.nan(dfs_1[[stat_list[i]]]) & !is.nan(dfs_2[[stat_list[i]]]))
  }

  
  # Compare magnitude and difference with each statistic between exp 1 and exp 2
  # Use decision rules that are defined in row_stats_toolbox, stored as attributes 
  # to dfs_1 and dfs_2
  for (i in seq_along(stat_list)) {
  
    # Determine if exp 1 has high null strength than exp 2
    # Candidates that returns NaNs are automatically designated as a wrong designation
    dfc[[paste("fract_", stat_list[i], "_1",strength_type, "2", sep='')]] <- 
      sum(match.fun(attr(dfs_1, strength_type)[[stat_list[i]]])
          (abs(dfs_1[[stat_list[i]]])[pseudo_indicies[,1]], 
            abs(dfs_2[[stat_list[i]]])[pseudo_indicies[,2]]), na.rm = TRUE) / df$n_samples
    # Calculate different of effect size between experiments
    dfc[[paste("mean_diff_", stat_list[i], "_1",strength_type,"2", sep='')]] <-
      abs(dfc[[paste("exp2_mean_", stat_list[i],sep='')]]) -
      abs(dfc[[paste("exp1_mean_", stat_list[i],sep='')]])
  }


  df_out <- cbind(df, dfc)
  
  if (any(is.nan(data.matrix(df_out))) || is.nan(dfc[[paste("fract_", stat_list[i], "_1",strength_type, "2", sep='')]]))
    { save(list = ls(all.names = TRUE),
    file = "temp/quantify_population_config.RData",envir = environment())}
  
  # Carry over data frame metadata
  user_attr <- attr(dfs_1,"user_attributes")
  for (i in seq_along(user_attr)) {attr(df_out,user_attr[i])<-attr(dfs_1,user_attr[i])}
  
  # Record average t-ratio for the samples generated
  # Simulations need to cover the parameter space
  # We define the parameter space is how "null" the results are, and use the
  # t-ratio to measure
  t_scores_1 <- row_tscore_2s(x_1a, x_1b)
  t_crits_1 <- qt(p = df$alpha_1, df = df$n_1a + df$n_1b -1, lower.tail = TRUE) 
  df_out$exp1_t_ratio_mean <- mean(t_scores_1/t_crits_1)
  
  if (!is.na(df$mu_2a)) {
  t_scores_2 <- row_tscore_2s(x_2a, x_2b)
  t_crits_2 <- qt(p = df$alpha_2, df = df$n_2a + df$n_2b - 1, lower.tail = TRUE) 
  df_out$exp2_t_ratio_mean <- mean(t_scores_2/t_crits_2)
  }
  
  # browser()
  # save(list = ls(all.names = TRUE), file = "temp/quantify_population_config.RData",envir = environment())
  # # load(file = "temp/quantify_population_config.RData")
  
  return(df_out)
}


tidy_strength_contest <- function (df, gt_colname, var_prefix, long_format = TRUE,
                                   ref_colname = NULL) {
  #' @description Converts comparison error rates for all candidate statistics 
  #' from wide format to long format. 
  #' Also applies normalization to error rates:
  #' 1) Subtract comparison error rate of candidate statistics from the ground 
  #' truth "error rate" determined by the underlying population parameters.
  #' 2) Normalize error rates from ground truth as a fraction of the error rate 
  #' from a priori selected gold standard reference statistic.
  #' 
  #' @param df input data
  #' @param gt_colname column names in df that serves as independent variable
  #'  and ground truth for determining whether exp 1 or 2 has higher 
  #'  strength for each pop. param set.
  #' @param var_prefix
  #' @param long_format flag to convert comparison error to long format
  #' @param ref_colname string of column name for candidate statistics that is a
  #'  gold standard reference variable that represent a best case comparison error 
  #'  rate. Error rates are normalized as a fraction of this best case.
  #' 
  #' @return df_tidy tidy (long) version of input dataframe
  
  save(list = ls(all.names = TRUE), file = "temp/tidy_strength_contest.RData",envir = environment())
  # # load(file = "temp/tidy_strength_contest.RData")
  
  # Check if gt_colname exists
  if (length(grep(gt_colname, names(df))) == 0) stop("gt_colname does not exist in dataframe")
  
  # Get list of all variables that match the variable suffix string
  matched_vars <- grep(var_prefix,colnames(df), perl=TRUE, value=TRUE)
  
  # Initialize new zeroed df with only selected vars
  data = matrix(0, dim(df)[1], length(matched_vars))
  colnames(data) <- matched_vars
  df_gt_sub <- as_tibble(data)
  
  # Fill in resulting values depending on variable type
  # Frac: fraction of samples for a given metric that satisfies a condition: to process,
  #        subtract these values from ground truth.
  # mean_diff: mean value of difference between groups:  to process,
  #        test for equivalent logical values
  for (n in seq(1,length(matched_vars), by = 1)) {
    if (var_prefix == "fract") {
      # If var_prefix is "fract.*", for cases where GT vector is TRUE, compute (1-x)
      df_gt_sub[[matched_vars[n]]][ df[[gt_colname]]] <- 
        1 - df[[matched_vars[n]]][ df[[gt_colname]]]
      df_gt_sub[[matched_vars[n]]][!df[[gt_colname]]] <-   
        df[[matched_vars[n]]][!df[[gt_colname]]]
    } else { # Matches mean_diff version of candidate statistics
      # if var_prefix is "mean_diff.*", subtract from GT
      df_gt_sub[matched_vars[n]] = df[matched_vars[n]] - gt_vector
    }
  }
  
  # If reference value included, subtract in paired fashion (reference value 
  # must be included in matched_variable list, with ground truth subtracted from it already)
  df_ref_sub <- df_gt_sub
  if (!is.null(ref_colname)) {stop("Reference normalization requested in tidy, discontinued")}
  #   for (n in seq(1,length(matched_vars), by = 1)) {
  #     df_ref_sub[matched_vars[n]] <- df_gt_sub[matched_vars[n]] -
  #       df_gt_sub[ref_colname]
  #   }
  #   print("Normalized by reference variable")
  # }
  
  # Flatten df_ref_sub into df_tidy (wide format to long format data frame)
  if (long_format) {
    df_tidy <- df_ref_sub %>% gather(matched_vars,factor_key=TRUE)
    names(df_tidy)[names(df_tidy) == "matched_vars"] <- "name"
  } else df_tidy <- df_ref_sub
  
  # Carry forward any user defined attributes
  user_attr <- attr(df,"user_attributes")
  for (i in seq_along(user_attr)) {attr(df_tidy,user_attr[i])<-attr(df,user_attr[i])}
  
  return(df_tidy)
}

pretty_stat_labels<- function(df, var_prefix) {
  #' @description Convert long levels names to shorter ones for pretty print in plots
  #' 
  #' @param df input dataframe
  #' @param var_prefix string prefix that matches the comparison error rate for 
  #' the candidate statistics
  #' 
  #' @return df_pretty converts levels of df into pretty format (for plotting)
  
  # save(list = ls(all.names = TRUE), file = "temp/pretty_stat_labels.RData",envir = environment())
  # load(file = "temp/pretty_stat_labels.RData")

  # Get current names for levels and the basenames for them to be matched against
  orig_levels <- levels(df$name)
  
  base_names <- paste(var_prefix,"_",attr(df,'varnames'),sep='')
  
  # Find index for each basename
  ind <- rep(0, length(base_names))
  for (n in seq(1,length(base_names))) {
    # new_levels_ind =1
    ind[n] <- which(regexpr(base_names[n],orig_levels)>0)
  }
  # Return new level names
  new_levels <- attr(df,"varnames_pretty")[ind]
  df_new_levels <- df
  levels(df_new_levels$name) <- unlist(unname(new_levels))

  return(df_new_levels)
}


plot_comparison_error <- function(df_pretty, fig_name, fig_path, y_ax_str) {
  #' @description Plot comparison error rates for candidates effect size statistics
  #' 
  #' @param df_pretty dataframe of comparison error rates in long pretty format
  #' @param fig_name filename of figure exported to disk
  #' @param fig_path path to output figures saved to disk
  #' @param y_ax_str String denoting independent variable for comparison error plot. 
  #' Either mu[DM], sigma[D], or df[D].
  #' 
  #' @return no return, exports figures to disk
  
  save(list = ls(all.names = TRUE), file = "temp/plot_comparison_error.RData",envir = environment())
  # load(file = "temp/plot_comparison_error.RData")
  
  # Calculate group means and corrected confidence intervals
  # Note: if it errors here with df_result having one group then plyr package 
  # was loaded before dplyr
  # browser();
  
  df_result <- df_pretty %>%   
    group_by(name) %>% 
    summarize(mean = mean(value), bs_ci_mean_str = toString(
      boot.ci(boot(value, function(x, ind)  mean(x[ind]), R = 10000),
              conf = 1-(0.05/choose(length(levels(df_pretty$name)), 2)),
              type = "basic" )$
        basic[c(4,5)]))
  df_result$bs_ci_mean_lower <- sapply(strsplit(df_result$bs_ci_mean_str,","), 
                                       function(x) as.numeric(x[1]))
  df_result$bs_ci_mean_upper <- sapply(strsplit(df_result$bs_ci_mean_str,","), 
                                       function(x) as.numeric(x[2]))
  # Check that group by name succeeded
  if (dim(df_result)[1]==1) {stop('dplr grouping failed, plyr loaded before dplyr?') }
 
  # If groups have no variance in mean, then boot returns NaNs, replace with mean
  df_result$bs_ci_mean_lower[is.na(df_result$bs_ci_mean_lower)] <- 
    df_result$mean[is.na(df_result$bs_ci_mean_lower)]
  df_result$bs_ci_mean_upper[is.na(df_result$bs_ci_mean_upper)] <- 
    df_result$mean[is.na(df_result$bs_ci_mean_upper)]
  
  df_result$is_mean_0.5 <- 0.5 >= df_result$bs_ci_mean_lower &
    0.5 <= df_result$bs_ci_mean_upper
  ci_range <-  c(min(df_result$bs_ci_mean_lower), max(df_result$bs_ci_mean_upper))
  
  # Export CSV table of means and relative means for results write-up
  mean_by_group <- df_result$mean
  mean_row <- matrix(rep(df_result$mean,length(levels(df_pretty$name))), 
                     ncol = length(levels(df_pretty$name)), 
                     nrow = length(levels(df_pretty$name)))
  rownames(mean_row) <- as.list(levels(df_pretty$name))
  mean_col <- t(mean_row)
  rmean_by_group <- (( mean_row - mean_col)/mean_col) * 100
  colnames(rmean_by_group) <- as.list(levels(df_pretty$name))
  names(mean_by_group) <- as.list(levels(df_pretty$name))
  # Calculate multiple comparison p value
  paired_results <- pairwise.t.test(df_pretty$value, df_pretty$name, p.adjust.method = "bonferroni",
                                    pool.sd = FALSE, paired = FALSE, alternative = c("two.sided"))
  paired_pvalues <- paired_results[[3]]
  # Write tables to file
  # Suppress warning with writing multiple tables to csv with column names
  csv_path <- paste(fig_path, "output_",str_replace(fig_name, ".tiff$", ".csv"), sep="")
  cat("Table 1: mean error rate for each group\n", file = csv_path)
  suppressWarnings(write.table(t(mean_by_group), csv_path, append = TRUE, 
                               col.names=TRUE, sep=","))
  cat("\nTable 2: relative mean error rate for each group\n",  append = TRUE, 
      file = csv_path)
  suppressWarnings(write.table(rmean_by_group, csv_path, col.names = TRUE, 
                               sep=",", append=TRUE))
  cat("\nTable 3: p value of mean error rate for each group\n",  append = TRUE, file = csv_path)
  suppressWarnings(write.table(paired_pvalues, csv_path, col.names = TRUE, 
                               sep=",", append=TRUE))
  
  # Labels of statistical significance for each group
  sig_labels = rep("",length(levels(df_pretty$name)))
  sig_colors = rep("black",length(levels(df_pretty$name)))
  sig_sizes = rep(4,length(levels(df_pretty$name)))
  siff_vjust = rep(0,length(levels(df_pretty$name)))
  # Set less than random to blue and -
  sig_labels[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  "-"
  sig_colors[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  
    rgb(47, 117, 181,maxColorValue = 255)
  sig_sizes[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  5
  siff_vjust[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  .017
  # Set greater than random to red and +
  sig_labels[df_result$bs_ci_mean_lower>0.5 & df_result$bs_ci_mean_upper>0.5] =  "+"
  sig_colors[df_result$bs_ci_mean_lower>0.5 & df_result$bs_ci_mean_upper>0.5] = 
    rgb(255, 0, 0,maxColorValue = 255)

  
  # Basic violin plot
  p <- ggplot(df_result, aes(x=name,  y=mean, group=name)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey") +
    geom_linerange(aes(ymin = bs_ci_mean_lower, ymax = bs_ci_mean_upper), size = 0.5) +
    geom_point(size=1,fill="white", shape = 1) + 
    xlab("Statistic") +
    ylab(parse(text= y_ax_str)) +
    scale_x_discrete(labels = parse(text = levels(df_pretty$name))) +
    expand_limits(y = c(0,1.1)) +
    geom_text(y = 1.07+siff_vjust, aes(label = sig_labels), 
              color = sig_colors, size = sig_sizes, vjust=0.5, hjust=0.5) +
    theme_classic() +  theme(text = element_text(size = 8))+
  scale_y_continuous(expand = c(0, 0))
  print(p)
  save_plot(paste(fig_path,  '/', fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 3, dpi = 600)

  return(df_result)
  
}


plot_mean_t_ratio_samples <- function(df_init, fig_path,fig_name) {
  
  save(list = ls(all.names = TRUE), file = "temp/plot_mean_t_ratio_samples.RData",envir = environment())
  # load(file = "temp/plot_mean_t_ratio_samples.RData")
  
  # Plot histogram of t_stat/t_critical to demonstrate how far from zero D is
  df <-tibble(group = as.factor(c(rep(1,dim(df_init)[1]),rep(2,dim(df_init)[1]))),
              t_ratio = c(df_init$exp1_t_ratio_mean, df_init$exp2_t_ratio_mean))
  p <- ggplot(data=df, aes(x = t_ratio, y = t_ratio, fill = group)) +
    geom_histogram(aes(y=stat(count / sum(count))), position = "identity", 
                   alpha=0.25, bins = 30) +
    geom_vline(xintercept = -1,linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    xlab( expression(bar(t)[stat]*phantom(.)/phantom(.)*abs(phantom(.)*t[crit]*phantom(.)))) +
    ylab( "Freq.") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01),expand = c(0, 0))
  ymax <- max(ggplot_build(p)$data[[1]]$ymax)
  # KS Test between both groups as they are plotted to see if hist. are different
  fill_id <- levels(as.factor(ggplot_build(p)$data[[1]]$fill))
  ks_p_val <- ks.test(subset(ggplot_build(p)$data[[1]], fill == fill_id[1])$y, 
                      subset(ggplot_build(p)$data[[1]], fill == fill_id[2])$y, 
                      alternative = "two.sided", exact = FALSE)
  p <- p + expand_limits(y = c(0, 1.1* ymax)) + 
    annotate("label",label=ifelse(ks_p_val$p.value>0.05, 
                                  sprintf('p = %.3f', ks_p_val$p.value),
                                  sprintf('p = %.2e', ks_p_val$p.value))
             , x=mean(range(df$t_ratio)), y=1.07*ymax, 
             size=2, fill = "white",label.size = NA)
  # print(p)
  save_plot(paste(fig_path, fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 1.2, dpi = 600)
  
}

process_strength_contest <- function(df_init, gt_colname, measure_pretty_str, out_path = paste(fig_path, "/temp",sep=''),
                                      fig_name, fig_path, var_prefix = "fract",
                                      parallel_sims = TRUE, is_plotted = TRUE, 
                                      stat_exclude_list= c("ldm", "rldm"),
                                      strength_type = "hnst", delta, is_delta_relative,
                                     use_pseudo_samples = FALSE) {
  #' @description 
  #' 
  #' @param df_init dataframe of population param sets by row
  #' @param gt_colname column names in df that serves as independent variable
  #'  and groundtruth for determining whether exp 1 or 2 has higher 
  #'  strength for each pop. param set.
  #' @param out_path path to export data files
  #' @param measure_pretty_str pretty format of independent variable
  #' @param fig_path path to export figures
  #' @param fig_name filename of figure exported to disk
  #' @param fig_path path to output figures saved to disk
  #' @param var_prefix string suffix to identify the columns that should use
  #' @param parallel_sims flat to process simulation (pop param sets) in parallel
  #' @param is_plotted flag whether results should be plotted and exported to disk
  #' @param stat_exclude_list list of candidate statistics to excluide in the 
  #' analysis
  #' 
  #' @return all_dfs a named list of all dataframes for each of the processing steps
  
  save(list = ls(all.names = TRUE), file = "temp/process_strength_contest.RData",envir = environment())
  # load(file = "temp/process_strength_contest.RData")
  dir.create(file.path(getwd(),out_path), showWarnings = FALSE)
  dir.create(file.path(getwd(),fig_path), showWarnings = FALSE)
  
  # Display ground truth fraction of E2>E1
  print(sprintf("%s (TRUE): %i", gt_colname, sum(df_init[[gt_colname]])))
  
  # Quantify effect sizes in untidy matrix
  df_es <- quantify_population_configs(df = df_init,overwrite = TRUE, out_path = out_path,
                                      data_file_name = paste(fig_name,".rds",sep = ""),
                                      parallel_sims = parallel_sims,
                                      stat_exclude_list=stat_exclude_list,
                                      strength_type = strength_type, delta = delta,
                                      is_delta_relative = is_delta_relative,
                                      use_pseudo_samples = use_pseudo_samples)
  
  is_single_exp = any(is.na(df_es$mu_2a)) || any(is.na(df_es$mu_2b))
  
  if (!is_single_exp) {
  
  # Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
  df_tidy <- tidy_strength_contest(df = df_es, gt_colname = gt_colname,
                                    var_prefix = var_prefix,long_format = TRUE,
                                    ref_colname = NULL)
  df_pretty <- pretty_stat_labels(df = df_tidy, var_prefix = var_prefix)
  
  # Plot effect size results
  if (is_plotted) {
    # df_compare_string <- tibble(lt="Lesser", gt = "Greater")
 

    df_inequality <- tibble("gt"="Greater", "lt"="Lesser")
    inequality_str <- df_inequality[[attr(df_init,paste('df_',strength_type,sep=''))[[gt_colname]]]]
    
    y_ax_str <- paste("Error~Rate~(~",inequality_str,"~phantom(.)*", measure_pretty_str,
          "*phantom(.))~phantom(.)~phantom(.)~phantom(.)~phantom(.)")
    
    df_plotted <- 
      plot_comparison_error(df = df_pretty, fig_name = fig_name, fig_path = fig_path,
                             y_ax_str = y_ax_str)
    
    plot_mean_t_ratio_samples(df = df_es, fig_path = fig_path, fig_name = paste("t-ratio_", fig_name, sep=""))
      
  }
  } else {df_tidy <- NULL; df_pretty = NULL; df_plotted = NULL;
  }
  
  # Package dataframes throughout processing into single list for return
  all_dfs <- vector(mode = "list", length = 4)
  names(all_dfs) <- c("df_es", "df_tidy", "df_pretty", "df_plotted")
  all_dfs[[1]] <- df_es; 
  all_dfs[[2]] <- df_tidy; 
  all_dfs[[3]] <- df_pretty; 
  if (exists("df_plotted")) {all_dfs[[4]] <- df_plotted; }
  
  # save(list = ls(all.names = TRUE), file = "temp/process_strength_contest.RData",envir = environment())
  
  # Some candidate stats may return NaNs, visualize wist histogram of fraction 
  # of NaN samples across simulations
  plot_nincluded_samples(df = df_es, var_prefix = "nincluded",
                               fig_name = fig_name, fig_path = fig_path)
  
  return(all_dfs)
}


plot_nincluded_samples<- function(df = df_es, var_prefix = "nincluded",fig_name = fig_name, fig_path = fig_path) {
  
  save(list = ls(all.names = TRUE), file = "temp/plot_nincluded_sampless.RData",envir = environment())
  # load(file = "temp/plot_nincluded_sampless.RData")
  
  
  # Get list of all variables that match the variable suffix string
  matched_vars <- grep(var_prefix,colnames(df), perl=TRUE, value=TRUE)
  
  
  df_nins <- df %>% pivot_longer(matched_vars, names_to = "variable", values_to = "count") %>% select(variable, count, n_actual_samples)
  df_nins$count_prc = 1-df_nins$count/df_nins$n_actual_samples
  # df_nins$variable <- sapply(df_nins$variable, function(x) {substring(x,17, nchar(x))})
  df_nins$variable <- as.factor(df_nins$variable)
  levels(df_nins$variable) <- as.factor(attr(df,"varnames_pretty") )
  
  gg <- ggplot(data = df_nins, aes(x = count_prc), fill = variable) +
    geom_histogram(bins = 20) + facet_wrap(~variable, labeller=label_parsed) +
    geom_vline(xintercept = 0.05, color = "red") +
    scale_y_continuous(trans='log10') +
    ylab("Log 10 Count") + xlab("fract discarded samples across sims") +
    theme_classic()
  print(gg)
  save_plot(paste(fig_path,  '/',   'histo_fract_discard', fig_name,
                  sep = ""), gg, ncol = 1, nrow = 1, 
            base_height = 5, base_width = 5, dpi = 600)

  # Check that all 
  if (any(df_nins$count_prc > 0.05))
  {save(list = ls(all.names = TRUE), file = "temp/plot_nincluded_sampless.RData",envir = environment());
    stop(sprintf('More than 5%% of samples in %i simulations returned NaN',
                  sum(df_nins$count_prc > 0.05)))}
  
  attr(df,"varnames_pretty")
}

plot_stats_covary_indvar <- function(df, indvar, indvar_pretty, fig_name, fig_path,
                                     dir_to_stronger = 1, alpha = 0.05) {
  #' @description Plot mean values across samples of candidate statistics versus 
  #' a swept independent variable, and calculates correlation with spearman rho.
  #' 
  #' @param df input dataframe of processed pop. params (returned from 
  #' process_strength_contest)
  #' @param indvar string of the column name describing the independent variable, 
  #' will be one of the strength measures
  #' @param fig_name base name of exported figure
  #' @param fig_path path of eported figure
  #' @param dir_to_stronger which direction leads to higher strength (+1,-1)
  #' @param alpha confidence level
  #' 
  #' @return df_mean_stat dataframe of spearman correlation and regression of each 
  #' candidate statistics versus the indepdent variable.
  
  
  save(list = ls(all.names = TRUE), file = "temp/plot_stats_covary_indvar.RData",envir = environment())
  # load(file = "temp/plot_stats_covary_indvar.RData")
  # browser();
  
  
  # Generate the names of all columns that refer to the mean and std of candidate statistics
  # Only candidate statistics from experiment 1 are used for this analysis, exp 2
  # are just dummy variables.
  exp1_mean_vars = paste("exp1_mean_", attr(df,"varnames"),sep="")
  exp1_sd_vars   = paste("exp1_sd_", attr(df,"varnames"),sep="")
  
  
  
  
  # Prepping data frame to plot how each the mean of each candidate statistic 
  # varies versus independent variable
  # __________________________________________________________________________
  # Calculate mean value of statistics for each value of indvar
  df$series <- seq(1, dim(df)[1],1)
  df_means <- df %>% gather("variable", "value", exp1_mean_vars)

  df_means <- pivot_longer(df,exp1_mean_vars, names_to = "variable", values_to = "mean_value") %>%
    select("series", c(all_of(indvar),"variable", "mean_value"))
  df_sd <- pivot_longer(df,exp1_mean_vars, names_to = "variable", values_to = "sd_value")
  # subset(df_means, variable = "exp1_mean_xdbar")
  df_means$variable <- factor(df_means$variable, levels = exp1_mean_vars)
  df_means$label <-  factor(as.character(attr(df,"varnames_pretty"))[match(df_means$variable,exp1_mean_vars)],
                            levels = attr(df,"varnames_pretty"))
  
  # The facet plot below would eror if the candidate statistics vary too much, 
  # so use log transform when needed
  # Log rescale values if any values are outside of specified range
  apply_log <- function(x) {if ((any(x<1e-5) || any(x>1e5)) && all(x>0)) {x=log10(x)} else{x=x}; return(x)}
  df_means2 = df_means %>% group_by(variable) %>% mutate(scale_mean_value = apply_log(mean_value))
  df_means2$is_scaled <- df_means2$mean_value != df_means2$scale_mean_value
  # Add log[10] to label for that group if scaling is applied
  df_is_scaled <- df_means2 %>% group_by(label) %>% summarize(is_scaled = all(is_scaled))
  levels(df_means2$label ) <-  as.character(ifelse(df_is_scaled$is_scaled, 
                                      paste("log[10](", attr(df,"varnames_pretty"),")",sep=""),
                                      attr(df,"varnames_pretty")))

  
  # Plot facet of how each candidate statistic changes
  gg <- ggplot(data = df_means2, aes_string(x = indvar, y = "scale_mean_value")) + 
    geom_point(size=.25) + facet_wrap(vars(label),labeller = label_parsed, ncol = 6, scales = "free_y") +
    ylab("Mean Value") + xlab(parse(text=indvar_pretty)) +
    # scale_y_continuous(breaks = NULL)+
    theme_classic(base_size=8) + 
    theme(strip.text.x = element_text( margin = margin( b = 1, t = 1) ))
  if (dir_to_stronger==1) {
    gg <- gg + scale_x_reverse(breaks=round(seq(min(df_means2[[indvar]]), max(df_means2[[indvar]]), 
     (max(df_means2[[indvar]]) - min(df_means2[[indvar]]))),2) )
  }else {
    gg<- gg+ scale_x_continuous(breaks=seq(min(df_means2[[indvar]]), max(df_means2[[indvar]]), 
                                           (max(df_means2[[indvar]]) - min(df_means2[[indvar]])))) 
  }
  print(gg)
  save_plot(paste(fig_path, '/', str_replace(fig_name,'\\..*$','_all_stats.png'), sep = ""), 
            gg, ncol = 1, nrow = 1, base_height = 1.75,  base_width = 6.25, dpi = 600)
 


  ## Examine correlation regression of candidate statistics to independent variable
  # Test: indvar must vary for correlation to be applied
  if (length(unique(df_means[[indvar]]))==1) {
    simpleError("Indepedent variable does not change, so cannot perform spearman correlation")
  }
  
  # Calculate spearman rho for the mean of each candidate stat versus independent var
  df_mean_stat = tibble(variable = exp1_mean_vars, 
                        spearman_rho=rep(NA,length(exp1_mean_vars)),
                        spearman_rho_low=rep(NA,length(exp1_mean_vars)), 
                        spearman_rho_high=rep(NA,length(exp1_mean_vars)), 
                        slope=rep(NA,length(exp1_mean_vars)), 
                        slope_low=rep(NA,length(exp1_mean_vars)),
                        slope_high=rep(NA,length(exp1_mean_vars)),
                        slope_nzero=rep(NA,length(exp1_mean_vars)))
  
  for (n in seq(1,length(exp1_mean_vars),1))  {
    # print(n)
    df_sub = subset(df_means, df_means$variable == exp1_mean_vars[n])
    
    #Calculate Spearman rho if stat varies, assign zero if not
    if (length(unique((abs(df_sub$mean_value)))) > 1) {
      # Spearman correlation
      tryCatch({
      ci_spearman = ci_cor(x = abs(df_sub[[indvar]])*-dir_to_stronger, 
                          y = abs(df_sub$mean_value), R = 1e5,
                          method = "spearman", type = "bootstrap", boot_type = "basic", 
                          probs = c(alpha/(2*length(exp1_mean_vars)), 
                                    1-alpha/(2*length(exp1_mean_vars))))
      ci_high = ci_spearman[[2]][2]
      ci_low = ci_spearman[[2]][1]
      ci_est = ci_spearman[[3]]
      }, error = function(e) {if (e$message=="length(ci) == 2L is not TRUE") 
        {ci_high <- 0; ci_low <- 0; ci_est <- 0}
        else {stop("plot_stats_covary_indvar:ci_cor: unrecognized error")
        }}
      )
      
      
      df_mean_stat$spearman_rho[n] = ci_est
      df_mean_stat$spearman_rho_low[n] = ci_low
      df_mean_stat$spearman_rho_high[n] = ci_high
    } else { 
      df_mean_stat$spearman_rho[n]    <- 0
      df_mean_stat$spearman_rho_low[n]  <- 0
      df_mean_stat$spearman_rho_high[n] <- 0
    }
    # Linear regression
    fit <- lm(y2 ~ y1, data = tibble(y1 = abs(df_sub[[indvar]])*-dir_to_stronger,
                                         y2 = abs(df_sub$mean_value)))
    conf_int <- confint(fit, "y1", level = 1-0.05/length(exp1_mean_vars))
    df_mean_stat$slope[n] = fit$coefficients[2]
    df_mean_stat$slope_low[n] = conf_int[1]
    df_mean_stat$slope_high[n] = conf_int[2]
  }
  # Sort variable names and add pretty column
  df_mean_stat <- df_mean_stat[match(exp1_mean_vars, df_mean_stat$variable),]
  df_mean_stat$label <- factor(attr(df,"varnames_pretty"), levels = attr(df,"varnames_pretty"))
  df_mean_stat$is_spearman_rho_sig <-  (df_mean_stat$spearman_rho_low<0 & df_mean_stat$spearman_rho_high<0) |
    (df_mean_stat$spearman_rho_low>0 & df_mean_stat$spearman_rho_high>0)
  df_mean_stat$is_slope_sig <-  (df_mean_stat$slope_low<0 & df_mean_stat$slope_high<0) |
    (df_mean_stat$slope_low>0 & df_mean_stat$slope_high>0)

  ## Plot 95% confidence interval of spearman correlation in annotated plot
  # __________________________________________________________________________
  # Plot notations for positive and negative spearman rho
  # Labels of statistical significance for each group
  sig_labels = rep("",length(attr(df,"varnames")))
  sig_colors = rep("black",length(attr(df,"varnames")))
  sig_sizes = rep(4,length(attr(df,"varnames")))
  siff_vjust = rep(0,length(attr(df,"varnames")))
  # Set less than random to blue and
  neg_ind <- sign(df_mean_stat$spearman_rho)==-1 & df_mean_stat$is_spearman_rho_sig
  pos_ind <- sign(df_mean_stat$spearman_rho)==1 & df_mean_stat$is_spearman_rho_sig
  
  sig_labels[neg_ind] =  "-"
  sig_colors[neg_ind] = rgb(47, 117, 181,maxColorValue = 255)
  sig_sizes[neg_ind] =  5
  siff_vjust[neg_ind] =  .017
  # Set greater than random to red and +
  sig_labels[pos_ind] =  "+"
  sig_colors[pos_ind] = rgb(255, 0, 0,maxColorValue = 255)
  
  # Plot spearman rho confidence intervals for each candidate stat  
  gg <- ggplot(data = df_mean_stat, aes(x = label, y = spearman_rho)) + 
    geom_point(shape=1, size=1) + ylim(-1,1) +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_linerange(aes(ymin = spearman_rho_low, ymax = spearman_rho_high)) +
    ylab(expression(Spearman*"'"*s~rho)) + xlab("Statistic") +
    scale_x_discrete(labels= parse(text = as.character(df_mean_stat$label))) +
    geom_text(y = 1.17+siff_vjust, aes(label = sig_labels), 
              color = sig_colors, size = sig_sizes, vjust=0.5, hjust=0.5) +
    coord_cartesian(ylim = c(-1,1.15), expand = TRUE) +
    theme_classic(base_size=8) + theme(legend.position="none") 
  print(gg)  
  save_plot(paste(fig_path, '/', fig_name, sep = ""), gg, ncol = 1, nrow = 1, 
            base_height = 1.75, base_asp = 3, base_width = 3, dpi = 600)
  
  
  # Plot values of of mu, rmu, sigma, rsigma of d and b over simulations, and df
  df_runs = tibble(Series = rep(seq(1,dim(df)[1],1),6),
                   param = c(rep("mu[DM]",dim(df)[1]), rep("r*mu[DM]",dim(df)[1]),
                             rep("sigma[D]",dim(df)[1]), rep("r*sigma[D]",dim(df)[1]),
                             rep("df[D]",dim(df)[1]), rep("alpha[DM]",dim(df)[1])), 
                   Value = c(df$mu_1dm, df$rmu_1dm, df$sigma_1d, df$rsigma_1d, 
                             df$df_1d, df$alpha_1)
  )
  df_runs$param <- factor(df_runs$param, levels = c("mu[DM]", "r*mu[DM]", "sigma[D]",
                                                    "r*sigma[D]","df[D]","alpha[DM]"))
  

  return(df_mean_stat)
}