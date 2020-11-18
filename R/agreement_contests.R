### Effect Size Contest



# TODO: for specifying offset distribution in generateData(), should the standard
# deviation for the control group be added by stds or variances?
# Load package manager
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
source("R/parallel_utils.R")
source("R/row_effect_sizes.R")
source("R/mmd.R")

# Parse all functions in file for parallel processing using user functions
row_effect_sizes_fun <- parse_functions_source("R/row_effect_sizes.R")
mmd_functions <- parse_functions_source("R/mmd.R")

# Dictionary to keep track of variables tracking effect size metrics
effect_size_dict <- vector(mode="list", length=4)
names(effect_size_dict) <- c("prefix", "base", "suffix","label")
effect_size_dict[[1]] <- c("fract", "mean_diff")
effect_size_dict[[2]] <- c("xdbar", "rxdbar", "sdmd", "rsdmd", "bf", "pvalue",
                           "tostp", "cohend", "mmd",  "rmmd","nrand")
effect_size_dict[[3]] <- c("d2gtd1","2m1")
effect_size_dict[[4]] <- c("bar(x)[DM]", "r*bar(x)[DM]", "s[DM]", "r*s[DM]","Bf", "p[NHST]*phantom(.)",
                           "~p[TOST]", "Cd" , "delta[M]",
                           "r*delta[M]","Rnd")


pop_params_from_aoffset <- function( n_samples, n_sims, 
                                    mus_1a, sigmas_1a, 
                                    mus_2a, sigmas_2a,
                                    mus_1ao, sigmas_1ao, 
                                    mus_2ao, sigmas_2ao,
                                    n_1a, n_1b, n_2a, n_2b) {
  # browser()
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
    sigma_2a = sigmas_2a, sigma_2b = sigmas_2b, sigma_2d = sigmas_2d
  )
  return(df)
}


pop_params_from_ab <- function( n_samples, n_sims, 
                                mus_1a, sigmas_1a, 
                                mus_2a, sigmas_2a,
                                mus_1b, sigmas_1b, 
                                mus_2b, sigmas_2b,
                                n_1a, n_1b, n_2a, n_2b) {
  # Calculate D based on A and B parameters
  mus_1d = mus_1b - mus_1a;  sigmas_1d = sqrt(sigmas_1a^2 + sigmas_1b^2)
  mus_2d = mus_2b - mus_2a;  sigmas_2d = sqrt(sigmas_2a^2 + sigmas_2b^2)
  # Initialize df with mu_1a, sigma_1a, mu_1b, sigma_1b, mu_1d,
  df = tibble( n_obs = n_obs, n_samples = n_samples,
               mu_1a = mus_1a, mu_1b = mus_1b, mu_1d = mus_1d,  n_1a = n_1a, n_1b = n_1b,
               mu_2a = mus_2a, mu_2b = mus_2b, mu_2d = mus_2d,  n_2a = n_2a, n_2b = n_2b, 
               sigma_1a = sigmas_1a, sigma_1b = sigmas_1b, sigma_1d = sigmas_1d,
               sigma_2a = sigmas_2a, sigma_2b = sigmas_2b, sigma_2d = sigmas_2d
  )
  return(df)
}


pop_params_switches <- function(df_init, switch_sign_mean_d, switch_sign_mean_ab, 
                                switch_group_ab, switch_exp_12) {
  df <- df_init
  # browser();
  # Randomly switch sign of D for both experiments, recalculate B
  if (switch_sign_mean_d) { mus_sign = sample(c(-1,1), n_sims, TRUE)
  df$mu_1d =  mus_sign * df$mu_1d
  df$mu_2d =  mus_sign * df$mu_2d
  
  df$mu_1b = df$mu_1a + df$mu_1d
  df$mu_2b = df$mu_2a + df$mu_2d
  }
  # Randomly switch sign of both group a and b for exp 1 and 2 separately, recalculate d
  if (switch_sign_mean_ab) {
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    df$mu_1a[switch_boolean] <- -df$mu_1a[switch_boolean]
    df$mu_2a[switch_boolean] <- -df$mu_2a[switch_boolean]
    df$mu_1b[switch_boolean] <- -df$mu_1b[switch_boolean]
    df$mu_2b[switch_boolean] <- -df$mu_2b[switch_boolean]
    
    # Recalculate D
    df$mu_1d = df$mu_1b - df$mu_1a; 
    df$mu_2d = df$mu_2b - df$mu_2a;
  }
  # Randomly switch assignment for A and B for EACH experiment, recalculate D
  if (switch_group_ab) {
    # Random switch binary vector
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    # Save temp variables for switch
    temp_mu_1a     <- df$mu_1a;    temp_sigma_1a  <- df$sigma_1a
    temp_mu_1b     <- df$mu_1b;    temp_sigma_1b  <- df$sigma_1b
    temp_mu_2a     <- df$mu_2a;    temp_sigma_2a  <- df$sigma_2a
    temp_mu_2b     <- df$mu_2b;    temp_sigma_2b  <- df$sigma_2b
    # Exp 1
    df$mu_1a[switch_boolean]      <- temp_mu_1b[switch_boolean]
    df$sigma_1a[switch_boolean]   <- temp_sigma_1b[switch_boolean]
    df$mu_1b[switch_boolean]      <- temp_mu_1a[switch_boolean]
    df$sigma_1b[switch_boolean]   <- temp_sigma_1a[switch_boolean]
    # Exp 2
    df$mu_2a[switch_boolean]      <- temp_mu_2b[switch_boolean]
    df$sigma_2a[switch_boolean]   <- temp_sigma_2b[switch_boolean]
    df$mu_2b[switch_boolean]      <- temp_mu_2a[switch_boolean]
    df$sigma_2b[switch_boolean]   <- temp_sigma_2a[switch_boolean]
    # Recalculate D
    df$mu_1d = df$mu_1b - df$mu_1a; df$sigma_1d = sqrt(df$sigma_1b^2 + df$sigma_1a^2)
    df$mu_2d = df$mu_2b - df$mu_2a; df$sigma_2d = sqrt(df$sigma_2b^2 + df$sigma_2a^2)
  }
  if (switch_exp_12) {
    # Save temp variables for switch
    temp_mu_1a     <- df$mu_1a;    temp_sigma_1a  <- df$sigma_1a
    temp_mu_1b     <- df$mu_1b;    temp_sigma_1b  <- df$sigma_1b
    temp_mu_2a     <- df$mu_2a;    temp_sigma_2a  <- df$sigma_2a
    temp_mu_2b     <- df$mu_2b;    temp_sigma_2b  <- df$sigma_2b
    # Determine which simulations to switch parameters for a and b
    switch_boolean <- sample(c(TRUE,FALSE), n_sims, TRUE)
    # Switch specified parameters
    df$mu_1a[switch_boolean]      <- temp_mu_2a[switch_boolean]
    df$sigma_1a[switch_boolean]   <- temp_sigma_2a[switch_boolean]
    df$mu_1b[switch_boolean]      <- temp_mu_2b[switch_boolean]
    df$sigma_1b[switch_boolean]   <- temp_sigma_2b[switch_boolean]
    
    df$mu_2a[switch_boolean]      <- temp_mu_1a[switch_boolean]
    df$sigma_2a[switch_boolean]   <- temp_sigma_1a[switch_boolean]
    df$mu_2b[switch_boolean]      <- temp_mu_1b[switch_boolean]
    df$sigma_2b[switch_boolean]   <- temp_sigma_1b[switch_boolean]
    # Recalculate D
    df$mu_1d = df$mu_1b - df$mu_1a; df$sigma_1d = sqrt(df$sigma_1b^2 + df$sigma_1a^2)
    df$mu_2d = df$mu_2b - df$mu_2a; df$sigma_2d = sqrt(df$sigma_2b^2 + df$sigma_2a^2)
  }
  return(df)
}


# default distribution for population parameters for Exp 1 {a,b}, Exp 2 {a,b}
generateExperiment_Data <- function(n_samples, n_sims, rand.seed,
                                    # Control group pop. parameters
                                    mus_1a, sigmas_1a, 
                                    mus_2a, sigmas_2a,
                                    # Experiment group pop. parameters
                                    mus_1b = NA, sigmas_1b = NA, 
                                    mus_2b = NA, sigmas_2b = NA,
                                    n_1a, n_1b, n_2a, n_2b,
                                    # Difference distribution pop. parameters
                                    mus_1ao = NA, sigmas_1ao = NA, 
                                    mus_2ao = NA, sigmas_2ao = NA,
                                    switch_group_ab = FALSE,
                                    switch_sign_mean_ab = FALSE,
                                    switch_sign_mean_d = FALSE,
                                    switch_exp_12 = FALSE,
                                    fig_name = "test.tiff",
                                    fig_path = "Figure/",
                                    gt_colnames, is_plotted = TRUE) {
  #' Generate simulated experiment data for two experiments 
  #' 
  #' @description Generate simulated experiment data for two experiments with 
  #' distributions for mu and sigma specified as functions
  #' 
  #' @param n_samples number of samples (collection of observations), number of 
  #' times that a simulated experiment is repeated
  #' @param n_obs number of measurements in a sample/experiment
  #' @param n_sims number of simulations run, which are sets of experiments with
  #'  a single set of parameters for each (mu, sigma, n_obs)
  #' @param rand.seed seed number of generating consistent random numbers
  #' @param mus1 vector specifying distribution for population mean for
  #'  experiment 1
  #' @param sigmas1 vector specifying distribution for population standard deviation for
  #'  experiment 1
  #' @param mus2 vector specifying distribution for population mean for
  #'  experiment 2
  #' @param sigmas2 vector specifying distribution for population mean for
  #'  experiment 12
  
  # browser();
  # Expand any singleton pop param arguments replicate to number of simulations
  input_args <- formalArgs(generateExperiment_Data)
  pargs <-grep("^(mus)|(sigmas)|(n_\\d)", input_args, value=TRUE)
  # For any pop param not equal in length to n_sims, expand
  for (n in seq_along(pargs)) {
    if (length(get(pargs[n]))==1) assign(pargs[n], rep(get(pargs[n]),n_sims))
  }

  # Record some parameter values for simulations
  set.seed(rand.seed +1)
  
  # Generate initial dataframe from params, no switching done yet  
  if (any(is.na(mus_1b))) {
    df_init <- 
      pop_params_from_aoffset( n_samples = n_samples, n_sims= n_sims,
                               mus_1a = mus_1a, sigmas_1a = sigmas_1a,  
                               mus_2a = mus_2a, sigmas_2a = sigmas_2a,
                               mus_1ao = mus_1ao, sigmas_1ao = sigmas_1ao, 
                               mus_2ao = mus_2ao, sigmas_2ao = sigmas_2ao,
                               n_1a = n_1a, n_2a = n_2a, n_1b = n_1b, n_2b = n_2b) 
  } else {
    df_init   <- 
      pop_params_from_ab( n_samples = n_samples, n_sims = n_sims, 
                          mus_1a = mus_1a, sigmas_1a = sigmas_1a,  mus_2a = mus_2a, sigmas_2a = sigmas_2a,
                          mus_1b = mus_1b, sigmas_1b = sigmas_1b, mus_2b = mus_2b, sigmas_2b = sigmas_2b,
                          n_1a = n_1a, n_2a = n_2a, n_1b = n_1b, n_2b = n_2b)
  }
  # Switch params if needed flag is specified
  df <- pop_params_switches(df = df_init, switch_sign_mean_d = switch_sign_mean_d, 
                            switch_sign_mean_ab = switch_sign_mean_ab, 
                            switch_group_ab = switch_group_ab,
                            switch_exp_12 = switch_exp_12)
  
  # Pop
  # Mean of difference and DM
  df$mu_1d <- df$mu_1b - df$mu_1a
  df$mu_2d <- df$mu_2b - df$mu_2a
  # Is: Exp2 mu[d] > Exp1 mu[d]
  df$is_mud_2gt1 <-  abs(df$mu_2d) > abs(df$mu_1d)
  
  # STD of the difference
  df$sigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2)
  df$sigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) 
  df$is_sigmad_2gt1 <-  df$sigma_2d   > df$sigma_1d
  
  # Degrees of freedom of the difference and DM
  df$df_1d <- df$n_1a + df$n_1b - 2
  df$df_2d <- df$n_2a + df$n_2b - 2
  df$is_dfd_2lt1 <- df$df_2d < df$df_1d
  df$is_dfdm_2lt1 <- df$df_2d < df$df_1d
  
  # Pooled standard deviation
  df$sigma_1pool <- sqrt( ( (df$n_1a-1)*df$sigma_1a^2 + (df$n_1b-1)*df$sigma_1b^2) /
                            (df$n_1a-1 + df$n_1b -1 )  )
  df$sigma_2pool <- sqrt( ( (df$n_2a-1)*df$sigma_2a^2 + (df$n_2b-1)*df$sigma_2b^2) /
                            (df$n_2a-1 + df$n_2b-1 ))
  df$is_sigmapool_2gt1 <- df$sigma_2pool > df$sigma_1pool
  
  # Mean and std of difference in means (taken from mean of D since we had the option
  # to invert the sign for D earlier in code)
  df$mu_1dm <- df$mu_1d
  df$mu_2dm <- df$mu_2d
  df$is_mudm_2gt1 <-  abs(df$mu_2dm) > abs(df$mu_1dm)
  
  # STD of the difference in means
  df$sigma_1dm <- sqrt(df$sigma_1a^2/n_1a + df$sigma_1b^2/n_1b)
  df$sigma_2dm <- sqrt(df$sigma_2a^2/n_2a + df$sigma_2b^2/n_2b)
  df$is_sigmadm_2gt1 <-  df$sigma_2dm > df$sigma_1dm
  
  # Calculate ratio of sigma_md/mu_md to determine how close DM is close to zero,
  # determines whether results are in null region of critical region of t-test
  df$mu_ov_sigma_1dm <- df$mu_1dm / df$sigma_1dm
  df$mu_ov_sigma_2dm <- df$mu_2dm / df$sigma_2dm

  # Statistics of difference of means distribution 
  df$rmu_1dm <- df$mu_1dm / df$mu_1a
  df$rmu_2dm <- df$mu_2dm / df$mu_2a
  df$is_rmudm_2gt1 <-  abs(df$rmu_2dm) > abs(df$rmu_1dm)

  # Relative sigma of difference
  df$rsigma_1d <- df$sigma_1d / abs(df$mu_1a + df$mu_1d/2)
  df$rsigma_2d <- df$sigma_2d / abs(df$mu_2a + df$mu_2d/2)
  df$is_rsigmad_2gt1 <-  df$rsigma_2d > df$rsigma_1d
  
  # Relative Pooled Sigma
  df$rsigma_1pool <- df$sigma_1pool / abs(df$mu_1a + df$mu_1d/2)
  df$rsigma_2pool <- df$sigma_2pool / abs(df$mu_2a + df$mu_2d/2)
  df$is_rsigmapool_2gt1 <-  df$rsigma_2pool > df$rsigma_1pool
  
  # sigma of the difference of means distirbution
  df$rsigma_1dm <- df$sigma_1dm / abs(df$mu_1a + df$mu_1dm/2)
  df$rsigma_2dm <- df$sigma_2dm / abs(df$mu_2a + df$mu_2dm/2)
  df$is_rsigmadm_2gt1 <-  df$rsigma_2dm > df$rsigma_1dm
  
  
  
  # Diff:  pop_mean2 - pop_mean1
  df$mean_mud_d2md1 <- df$mu_2dm - df$mu_1dm
  df$mean_rmud_d2md1 <- (df$mu_2dm/df$mu_2a) - (df$mu_1dm/df$mu_1a)
  # Diff:  pop_std2 - pop_std1
  df$mean_sigmadm_2m1 <- df$sigma_2dm - df$sigma_1dm
  df$mean_rsigmadm_2m1 <- df$sigma_2dm/df$mu_2a - df$sigma_1dm/df$mu_1a
  
  # Append columns for effect sizes, since multiple columns are used to analyze
  # each effect size, a dictionary of prefix, base, and suffix variable names 
  # are used.
  df[ paste(effect_size_dict$prefix[1], effect_size_dict$base,
            effect_size_dict$suffix[1], sep="_") ] <- rep(NaN,n_sims)
  df[ paste(effect_size_dict$prefix[2], effect_size_dict$base,
            effect_size_dict$suffix[2], sep="_") ] <- rep(NaN,n_sims)
  
  # Plot generated population parameters
  if (is_plotted){
    plot_population_params(df, fig_name = fig_name, fig_path = fig_path, 
                           gt_colnames = gt_colnames)
  } else {
    # Plot values of of mu, rmu, sigma, rsigma of d and b over simulations, and df
    df_runs = tibble(Series = rep(seq(1,dim(df)[1],1),5),
                     param = c(rep("mu[DM]",dim(df)[1]), rep("r*mu[DM]",dim(df)[1]),
                               rep("sigma[pool]",dim(df)[1]), rep("r*sigma[pool]",dim(df)[1]),
                               rep("df[pool]",dim(df)[1])), 
                     Value = c(df$mu_1dm, df$rmu_1dm, df$sigma_1d, df$rsigma_1d, df$df_1d)
    )
    df_runs$param <- factor(df_runs$param, levels = c("mu[DM]", "r*mu[DM]", "sigma[pool]","r*sigma[pool]","df[pool]"))

    gg <- ggplot(data = df_runs, aes(x = Series, y = Value)) +
      geom_line() +
      facet_wrap(vars(param), nrow=3,ncol=2,scales="free_y",labeller=label_parsed) +
      theme_classic(base_size = 8) +
      theme(strip.text.x = element_text( margin = margin( b = 0, t = 0) ),
            panel.spacing = unit(0, "lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(fill = NA,colour = "black")) 
    print(gg)
    save_plot(paste(fig_path, str_replace(fig_name,"\\.[a-z]*$","_params.tiff"), sep = ""), gg, ncol = 1, nrow = 1, 
              base_height = 1.8, base_asp = 4, base_width = 3, dpi = 600)
  }
  return(df)
}


plot_population_params <- function(df_init, gt_colnames,fig_name,fig_path){
  #'
  #'
  #'
  #'
  #'
  #'
  #'
  #'
  # save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  # load(file = "temp/debug.RData")
  
  # Output csv of agreement of input parameters to each individual input parameter
  param_fields = c("is_mudm_2gt1","is_rmudm_2gt1","is_sigmad_2gt1",
                   "is_rsigmad_2gt1", "is_dfdm_2lt1")
  
  bv_gt_colnames <- sapply(gt_colnames, function(x) any(x==param_fields))
  if (!all(sapply(gt_colnames, function(x) any(x==param_fields)))) {stop("gt_colnames is not correctly worded")}
  
  
  # Calculate indices of colname
  gt_param_inds <- unname(sapply(gt_colnames,function(x){pmatch(x,param_fields)}))
  gt_param_labels <- c("abs(~mu[DM]*phantom(.))",  "abs(~r*mu[DM]*phantom(.))",
                       "~sigma[pool]*phantom(.)", "r*sigma[pool]", "df[pool]")

  # Calculate agreement matrix: test if nonrandom agreement between each parameter
  # versus every other.
  n_agreement = matrix(0, ncol = length(param_fields), nrow = length(param_fields))
  pwise_binom_p <- n_agreement
  
  for (r in seq(1,length(param_fields),1)) {
    for (c in seq(1,length(param_fields),1)) {
      n_agreement[r,c]   <- sum(df_init[[param_fields[r]]] == df_init[[param_fields[c]]])
      pwise_binom_p[r,c] <-
        prop.test(n_agreement[r,c], dim(df_init)[1], alternative = "two.sided",
                  conf.level = 1-0.05/(5*length(gt_colnames)), correct = TRUE)$p.value
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
  df_params <- tibble(group = c("mu", "rmu", "sigma", "rsigma", "df"), 
                      estimate = rep(0,length(param_fields)), lci = rep(0,length(param_fields)),
                      uci = rep(0,length(param_fields)))
  for (n in seq_along(param_fields)) {
    binom <- prop.test(sum(df_init[[param_fields[n]]]), dim(df_init)[1], 
                                  conf.level=1-0.05/(5*length(gt_colnames)), correct = FALSE)
    df_params$estimate[n]  <-binom$estimate
    df_params$lci[n]       <-binom$conf.int[1]
    df_params$uci[n]       <-binom$conf.int[2]
  }
  df_params$group <- factor(df_params$group, levels = c("mu", "rmu", "sigma", "rsigma", "df"))
  
  
  # Plot confidence interval of success rate for each parameter and their agreement
  gg <- ggplot(df_params, aes(x = group,  y = estimate)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey",linetype="dashed") +
    geom_linerange(aes(ymin = lci, ymax = uci), size = 0.5) +
    geom_point(size = 1.25, fill = "white", shape = 1) + 
    ylab(expression(Frac.~Exp[1]~hat.~Exp[2]~~~~~~~~~~~~~phantom("."))) +
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
    save_plot(paste(fig_path, 'gt_',fig_name, ".tiff", sep = ""), gg, ncol = 1, nrow = 1, 
              base_height = 1.5, base_asp = 3, base_width = 2, dpi = 600)
    
  }
  print(gg)
  save_plot(paste(fig_path, 'gt_',fig_name, ".tiff", sep = ""), gg, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 2, dpi = 600)
  

  # Export csv file for agreement between each variable to others
  # Plot histogram of mu[D]/sigma[D] to demonstrate how far from zero D is  
  df <-tibble(group = as.factor(c(rep(1,dim(df_init)[1]),rep(2,dim(df_init)[1]))),
              mu_ov_sigma = c(df_init$mu_1dm/df_init$sigma_1dm,
                              df_init$mu_2dm/df_init$sigma_2dm))
  p <- ggplot(df, aes(x = mu_ov_sigma, y = mu_ov_sigma, fill = group)) +
    geom_histogram(aes(y=stat(count / sum(count))), position = "identity", 
                   alpha=0.25, bins = 30) +
    geom_vline(xintercept = -2.6,linetype = "dashed") +
    geom_vline(xintercept = 2.6, linetype = "dashed") +
    xlab( expression(mu[DM]*phantom(.)/phantom(.)*sigma[DM])) +
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
               , x=0, y=1.07*ymax, 
             size=2, fill = "white",label.size = NA)
  # print(p)
  save_plot(paste(fig_path, 'mu_ov_sigma_',fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 1.2, dpi = 600)
  
}



quantify_esize_simulation <- function(df, include_bf = FALSE, rand.seed = 0, 
                                      parallelize_bf = FALSE) {
  if (dim(df)[1] != 1) stop("Need to input a single row of df")
  set.seed(rand.seed)

  # Use Exp 1 and 2 coefficients to generate data from normalized base data
  x_1a = matrix(rnorm(df$n_samples * df$n_1a, mean = df$mu_1a, 
                      sd = df$sigma_1a), nrow = df$n_samples, 
                ncol = df$n_1a)
  x_1b = matrix(rnorm(df$n_samples * df$n_1b, mean = df$mu_1b, 
                      sd = df$sigma_1b), nrow = df$n_samples, 
                ncol = df$n_1b)
  
  x_2a = matrix(rnorm(df$n_samples * df$n_2a, mean = df$mu_2a, 
                      sd = df$sigma_2a), nrow = df$n_samples, 
                ncol = df$n_2a)
  x_2b = matrix(rnorm(df$n_samples * df$n_2b, mean = df$mu_2b, 
                      sd = df$sigma_2b), nrow = df$n_samples, 
                ncol = df$n_2b)
  
  
  # Means
  xdbar_1d = rowMeans(x_1b) - rowMeans(x_1a)
  xdbar_2d = rowMeans(x_2b) - rowMeans(x_2a)
  df$exp1_mean_xdbar = mean(xdbar_1d)
  df$exp2_mean_xdbar =  mean(xdbar_2d)
  df$exp1_sd_xdbar = sd(xdbar_1d)
  df$exp2_sd_xdbar =  sd(xdbar_2d)
  diff_xdbar = abs(xdbar_2d) - abs(xdbar_1d)
  df$fract_xdbar_d2gtd1 = sum(diff_xdbar > 0) / df$n_samples
  df$mean_diff_xdbar_2m1  = mean(diff_xdbar)
  
  # Stds
  s_1md = sqrt(rowSds(x_1a)^2/df$n_1a + rowSds(x_1b)^2/df$n_1b)
  s_2md = sqrt(rowSds(x_2a)^2/df$n_2a + rowSds(x_2b)^2/df$n_2b)
  df$exp1_mean_sdmd = mean(s_1md)
  df$exp2_mean_sdmd = mean(s_2md)
  df$exp1_sd_sdmd   = sd(s_1md)
  df$exp2_sd_sdmd   = sd(s_2md)
  diff_sd = s_2md - s_1md
  df$fract_sdmd_d2gtd1  = sum(diff_sd > 0) / df$n_samples
  df$mean_diff_sdmd_2m1 = mean(diff_sd)
  
  # Rel Means: mean divided by control mean
  exp1_rxdbar = xdbar_1d/rowMeans(x_1a)
  exp2_rxdbar = xdbar_2d/rowMeans(x_2a)
  df$exp1_mean_rxdbar = mean(exp1_rxdbar)
  df$exp2_mean_rxdbar = mean(exp2_rxdbar)
  df$exp1_sd_rxdbar   = sd(exp1_rxdbar)
  df$exp2_sd_rxdbar   = sd(exp2_rxdbar)
  diff_rxdbar = abs(exp2_rxdbar) - abs(exp1_rxdbar)
  df$fract_rxdbar_d2gtd1 = sum( diff_rxdbar > 0) / df$n_samples
  df$mean_diff_rxdbar_2m1  =  mean(diff_rxdbar)
  
  # Rel STDs: sd divided by control mean
  # Relaltive STD is halfway between xbar_A and xbar_B
  exp1_rsd_md = s_1md / (rowMeans(x_1a) + 0.5 * xdbar_1d)
  exp2_rsd_md = s_2md / (rowMeans(x_2a) + 0.5 * xdbar_2d)
  df$exp1_mean_rsdmd = mean(exp1_rsd_md)
  df$exp2_mean_rsdmd = mean(exp2_rsd_md)
  df$exp1_sd_rsdmd   = sd(exp1_rsd_md)
  df$exp2_sd_rsdmd    = sd(exp2_rsd_md)
  diff_rsd = exp2_rsd_md - exp1_rsd_md
  df$fract_rsdmd_d2gtd1 = sum( diff_rsd > 0) / df$n_samples
  df$mean_diff_rsdmd_2m1  = mean(diff_rsd)
  
  # Bayes factor
  if (include_bf) {
    bf_1 = row_bayesf_2s(x_1a, x_1b, parallelize = parallelize_bf, paired = FALSE)
    bf_2 = row_bayesf_2s(x_2a, x_2b, parallelize = parallelize_bf, paired = FALSE)
    df$exp1_mean_bf = mean(bf_1)
    df$exp2_mean_bf = mean(bf_2)
    df$exp1_sd_bf   = sd(bf_1)
    df$exp2_sd_bf   = sd(bf_2)
    diff_abs_bf <- abs(bf_2) - abs(bf_1)
    df$fract_bf_d2gtd1 = sum(diff_abs_bf > 0) /
      df$n_samples
    df$mean_diff_bf_2m1 = mean(diff_abs_bf)
  } else {
    df$fract_bf_d2gtd1 = 0
    df$mean_diff_bf_2m1 = 0
  }
  
  # NHST P-value 
  # The more equal experiment will have a larger p-value
  z_score_1 <- row_zscore_2s(x_1b, x_1a)
  z_score_2 <- row_zscore_2s(x_2b, x_2a)
  pvalue_1 = 2*pnorm(-abs(z_score_1))
  pvalue_2 = 2*pnorm(-abs(z_score_2))
  df$exp1_mean_pvalue = mean(pvalue_1)
  df$exp2_mean_pvalue = mean(pvalue_2)
  df$exp1_sd_pvalue = sd(pvalue_1)
  df$exp2_sd_pvalue = sd(pvalue_2)
  diff_pvalue <-  pvalue_2 - pvalue_1
  # Must switch signs since larger is more similiar
  df$fract_pvalue_d2gtd1 = sum(-diff_pvalue > 0) / df$n_samples
  df$mean_diff_pvalue_2m1 = mean(diff_pvalue)
  
  # TOST p value (Two tailed equivalence test)
  tostp_1 <- row_tost_2s(x_1b, x_1a,low_eqbound = -.1,high_eqbound = .1)
  tostp_2 <- row_tost_2s(x_2b, x_2a,low_eqbound = -.1,high_eqbound = .1)
  df$exp1_mean_tostp = mean(tostp_1)
  df$exp2_mean_tostp = mean(tostp_2)
  df$exp1_sd_tostp = sd(tostp_1)
  df$exp2_sd_tostp = sd(tostp_2)
  diff_tostp <- tostp_2 - tostp_1
  df$fract_tostp_d2gtd1 <- sum(diff_tostp > 0) / df$n_samples
  df$mean_diff_tostp_2m1 <- mean(diff_tostp) 
  
  # Cohens D
  cohend_1 = row_cohend(x_1a, x_1b)
  cohend_2 = row_cohend(x_2a, x_2b) 
  df$exp1_mean_cohend = mean(cohend_1)
  df$exp2_mean_cohend = mean(cohend_2)
  df$exp1_sd_cohend = sd(cohend_1)
  df$exp2_sd_cohend = sd(cohend_2)
  diff_cohend = abs(cohend_2) - abs(cohend_1)
  df$fract_cohend_d2gtd1 = sum(diff_cohend > 0) / df$n_samples
  df$mean_diff_cohend_2m1 =  mean(diff_cohend)
  
  # Most Mean Diff
  mmd_1 = row_mmd_2s_zdist(x_1a, x_1b)
  mmd_2 = row_mmd_2s_zdist(x_2a, x_2b)
  df$exp1_mean_mmd = mean(mmd_1)
  df$exp2_mean_mmd = mean(mmd_2)
  df$exp1_sd_mmd = sd(mmd_1)
  df$exp2_sd_mmd = sd(mmd_2)
  diff_most_mean_diff = mmd_2 - mmd_1
  df$fract_mmd_d2gtd1 = sum(diff_most_mean_diff > 0) / df$n_samples
  df$mean_diff_mmd_2m1 = mean(diff_most_mean_diff)
  
  # Relative Most Mean Diff
  rmmd_1 = mmd_1 / rowMeans(x_1a)
  rmmd_2 = mmd_2 / rowMeans(x_2a)
  df$exp1_mean_rmmd = mean(rmmd_1)
  df$exp2_mean_rmmd = mean(rmmd_2)
  df$exp1_sd_rmmd = sd(rmmd_1)
  df$exp2_sd_rmmd = sd(rmmd_2)
  diff_rmmd =  rmmd_2 - rmmd_1
  df$fract_rmmd_d2gtd1 = sum(diff_rmmd > 0) / df$n_samples
  df$mean_diff_rmmd_2m1 = mean(diff_rmmd)
  
  
  # Random group
  nrand_1 = rowMeans(matrix(rnorm(df$n_samples * df$df_1d, mean = 0, sd = 1), 
                            nrow = df$n_samples, ncol = df$df_1d))
  nrand_2 = rowMeans(matrix(rnorm(df$n_samples * df$df_2d, mean = 0, sd = 1), 
                            nrow = df$n_samples, ncol = df$df_2d))
  df$exp1_mean_nrand = mean(nrand_1)
  df$exp2_mean_nrand = mean(nrand_2)
  df$exp1_sd_nrand = sd(nrand_1)
  df$exp2_sd_nrand = sd(nrand_2)
  diff_nrand = nrand_2 - nrand_1
  df$fract_nrand_d2gtd1 = sum(diff_nrand > 0 ) / df$n_samples
  df$mean_diff_nrand_2m1 = mean(diff_nrand)
    
  # browser()
  return(df)
}

  
  
quantify_esize_simulations <- function(df_in, overwrite = TRUE,
                                out_path = "temp/", data_file_name,
                                rand.seed = 0, include_bf = TRUE, 
                                parallel_sims = TRUE) {
  #' Simulate experiments generated from generateExperiment_Data() and calculates
  #'  various effect sizes
  #' 
  #' @description Given a data frame of parameters and input data for each 
  #' experiment group, quantifies data with various effect size metrics.
  #' 
  #' @param df input data frame generated from generateExperiment_Data() that 
  #' holds parameters for simulated data, along with initialized fields to store
  #'  effect size metrics
  #' @param x_a base input data from group a (gets transformed into actual input
  #'  data with parameters within df)
  #' @param x_b base input data from group b (gets transformed into actual input
  #'  data with parameters within df)
  #' @param out_path file path to save results to disk
  #' @param overwrite if results file already exists in out_path, skip 
  #' calculation and load from disk

  # Initialize output data frame
  # parallel_sims = FALSE
  n_sims = dim(df_in)[1]
  df <- df_in
  
  # browser();
  # save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  # # load(file = "temp/debug.RData")
  
  
  # Only perform simulations if results not saved to disk
  if (!file.exists(paste(out_path,data_file_name,sep="")) | overwrite) {
    
    if (parallel_sims) { print("starting parallel processing")
      # Setup parallel back end to use many processors
      cl = makeCluster(detectCores()[1]-1)
      registerDoParallel(cl)
      
      df <- foreach(n = 1:n_sims, .combine = rbind, .packages = 
                      c("BayesFactor","TOSTER"),
                .export = c(row_effect_sizes_fun, mmd_functions,
                            "quantify_esize_simulation")) %dopar% {
        #calling a function
        tempMatrix <- quantify_esize_simulation(df[n,], include_bf, 
                                               rand.seed = rand.seed+n,
                                               parallelize_bf = FALSE) 
        tempMatrix
      }
      #stop cluster
      stopCluster(cl)
    }else{
      df_list = list();
      # Process effect sizes serially and bind rows into one dataframe
      for (n in seq(1,n_sims,1)) {
        df_list[[n]] <- quantify_esize_simulation(df_in[n,], include_bf, rand.seed = rand.seed+n, 
                                            parallelize_bf = FALSE) 
      }
      df <- bind_rows(df_list)
    }
    
    # Save dataframe results to a file
    saveRDS(df, file = paste(out_path, data_file_name, sep=""))
  } else {
    # Restore the dataframe results from disk
    df <- readRDS(file = paste(out_path, data_file_name,sep=""))
  }
  
  return(df)
}



tidy_esize_simulations <- function (df, gt_colname, var_suffix, long_format = TRUE,
                                   ref_colname = NULL) {
  #' Normalize a subset of variables in df and then flatten data frame
  #' 
  #' @description Normalize a subset of variables in df, subract/compare to groundtruth, and 
  #' then subtract out the reference value if included
  #' 
  #' @param df
  #' @param gt_colname
  #' @param var_suffix
  #' @param long_format
  #' @param ref_colname
  
  # Check if gt_colname exists
  if (length(grep(gt_colname, names(df))) == 0) stop("gt_colname does not exist in dataframe")
  
  # Get list of all variables that match the variable suffix string
  matched_vars <- grep(var_suffix,colnames(df), perl=TRUE, value=TRUE)
  
  # Initialize new zeroed df with only selected vars
  data = matrix(0, dim(df)[1], length(matched_vars))
  colnames(data) <- matched_vars
  df_gt_sub <- as_tibble(data)
  
  # Fill in resulting values depending on variable type
  # Frac: fraction of metrics of samples that satisfy some condition: to process,
  #        subtract these values from ground truth.
  # mean_diff: mean value of difference between groups:  to process,
  #        test for equaivalent logical values
  for (n in seq(1,length(matched_vars), by = 1)) {
    if (var_suffix == "fract") {
      # If var_suffix is "fract.*", for cases where GT vector is TRUE, compute (1-x)
      df_gt_sub[[matched_vars[n]]][ df[[gt_colname]]] <- 
        1 - df[[matched_vars[n]]][ df[[gt_colname]]]
      df_gt_sub[[matched_vars[n]]][!df[[gt_colname]]] <-   
        df[[matched_vars[n]]][!df[[gt_colname]]]
    } else {
      # if var_suffix is "mean_diff.*", subtract from GT
      df_gt_sub[matched_vars[n]] = df[matched_vars[n]] - gt_vector
    }
  }
  
  # If reference value included, subtract in paired fashion (reference value 
  # must be included in matched_variable list, with groundtruth subtract from it)
  df_ref_sub <- df_gt_sub
  if (!is.null(ref_colname)) {
    for (n in seq(1,length(matched_vars), by = 1)) {
      df_ref_sub[matched_vars[n]] <- df_gt_sub[matched_vars[n]] - 
        df_gt_sub[ref_colname]
    }
  }
  
  # Flatten df_ref_sub into df_tidy (wide format to long format data frame)
  if (long_format) {
    df_tidy <- df_ref_sub %>% gather(matched_vars,factor_key=TRUE)
    names(df_tidy)[names(df_tidy) == "matched_vars"] <- "name"
  } else df_tidy <- df_ref_sub
  
  # browser()
  return(df_tidy)
}

pretty_esize_levels<- function(df,base_names, pretty_names, var_suffix) {
  #' Convert long levels names to shorter ones for pretty print in plots
  #' 
  #' @description 
  #' @param 
  #' @param 
  #' @param 
  #' 
  
  # Get current names for levels and the basenames for them to be matched against
  orig_levels <- levels(df$name)
  # base_names <- paste(var_suffix,"_", base_names, "_", sep="") 
  
  # Find index for each basename
  ind <- rep(0, length(base_names))
  for (n in seq(1,length(base_names))) {
    new_levels_ind =1
    ind[n] <- which(regexpr(base_names[n],orig_levels)>0)
  }
  # Return new level names
  new_levels <- effect_size_dict$label[ind]
  
  # df_new_levels <- df
  # levels(df_new_levels$name) <- new_levels
  levels(df$name) <- new_levels
  
  return(df)
  
}


plot_esize_simulations <- function(df_pretty, fig_name, fig_path, y_ax_str, comp_dir = "Lower") {
  
  # Calculate group means and corrected confidence intervals
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
  sig_colors[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  rgb(47, 117, 181,maxColorValue = 255)
  sig_sizes[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  5
  siff_vjust[df_result$bs_ci_mean_lower<0.5 & df_result$bs_ci_mean_upper<0.5] =  .017
  # Set greater than random to red and +
  sig_labels[df_result$bs_ci_mean_lower>0.5 & df_result$bs_ci_mean_upper>0.5] =  "+"
  sig_colors[df_result$bs_ci_mean_lower>0.5 & df_result$bs_ci_mean_upper>0.5] =  rgb(255, 0, 0,maxColorValue = 255)

  # browser()
  # Basic violin plot
  p <- ggplot(df_result, aes(x=name,  y=mean, group=name)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey") +
    geom_linerange(aes(ymin = bs_ci_mean_lower, ymax = bs_ci_mean_upper), size = 0.5) +
    geom_point(size=1,fill="white", shape = 1) + 
    xlab("Statistic") +
    ylab(parse(text=paste("Error~Rate~(~",comp_dir,"~phantom(.)*", y_ax_str,"*phantom(.))~phantom(.)~phantom(.)~phantom(.)~phantom(.)"))) +
    scale_x_discrete(labels = parse(text = levels(df_pretty$name))) +
    expand_limits(y = c(0,1.1)) +
    geom_text(y = 1.07+siff_vjust, aes(label = sig_labels), 
              color = sig_colors, size = sig_sizes, vjust=0.5, hjust=0.5) +
    theme_classic() +  theme(text = element_text(size = 8))+
  scale_y_continuous(expand = c(0, 0))
  print(p)
  save_plot(paste(fig_path, fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 3, dpi = 600)
   #browser()

  return(df_result)
  
}

process_esize_simulations <- function(df_init, gt_colname, y_ax_str, out_path = "temp/",
                                      fig_name, fig_path, var_suffix = "fract",include_bf = TRUE,
                                      parallel_sims = TRUE, is_plotted = TRUE, comp_dir = "Lower") {
  # browser();
  
  # Display ground truth fraction of E2>E1
  print(sprintf("%s (TRUE): %i", gt_colname, sum(df_init[[gt_colname]])))
  
  # Quantify effect sizes in untidy matrix
  df_es <- quantify_esize_simulations(df = df_init,overwrite = TRUE, out_path = out_path,
                                      data_file_name = paste(fig_name,".rds",sep = ""),
                                      include_bf = include_bf,parallel_sims = parallel_sims)
  # Tidy matrix by subtracting ground truth and normalizing to a reference variable if necessary
  df_tidy <- tidy_esize_simulations(df = df_es, gt_colname = gt_colname,
                                    var_suffix = var_suffix,long_format = TRUE,
                                    ref_colname = NULL)
  df_pretty <- 
    pretty_esize_levels(df = df_tidy, base_names = paste(var_suffix,"_",
                                                         effect_size_dict$base, "_", sep=""),
                        pretty_names = effect_size_dict$label, 
                        var_suffix = var_suffix)
  
  # Plot effect size results
  if (is_plotted) {
    df_plotted <- plot_esize_simulations(df = df_pretty, fig_name = fig_name, 
                                         fig_path = fig_path, y_ax_str = y_ax_str, comp_dir = comp_dir)
  }
  
  # Package dataframes throughout processing into single list for return
  all_dfs <- vector(mode = "list", length = 4)
  names(all_dfs) <- c("df_es", "df_tidy", "df_pretty", "df_plotted")
  all_dfs[[1]] <- df_es; 
  all_dfs[[2]] <- df_tidy; 
  all_dfs[[3]] <- df_pretty; 
  if (exists("df_plotted")) {all_dfs[[4]] <- df_plotted; }
  return(all_dfs)
}



lineplot_indvar_vs_stats <- function(df, indvar, fig_name, fig_path, stats_basenames = effect_size_dict[[2]], 
                                     stats_labels = effect_size_dict[[4]], dir_to_agreement=1, alpha = 0.05) {
  #' Plot correlation of independent variable versus mean values of all statistics
  #' 
  #' 
  #' 
  #' 
  save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  # load(file = "temp/debug.RData")
  
  # Filter for stats metrics
  col_list <- colnames(df)
  exp1_mean_vars = paste("exp1_mean_", stats_basenames,sep="")
  exp1_sd_vars   = paste("exp1_sd_", stats_basenames,sep="")
  
  # Calculate mean value of statistics for each value of indvar
  df_means <- df %>% gather("variable", "value", exp1_mean_vars)
 
  df_means <- pivot_longer(df,exp1_mean_vars, names_to = "variable", values_to = "mean_value") %>%
    select(c(all_of(indvar),"variable", "mean_value"))
  df_sd <- pivot_longer(df,exp1_mean_vars, names_to = "variable", values_to = "sd_value")
  # subset(df_means, variable = "exp1_mean_xdbar")
  df_means$variable <- as.factor(df_means$variable)
  
  
  if (length(unique(df_means[[indvar]]))==1) {
    simpleError("Indepedent variable does not change, so cannot perform pearson")
  }

  save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  # load(file = "temp/debug.RData")
  
  
  df_mean_stat = tibble(variable = exp1_mean_vars, pearson_rho=rep(NA,length(exp1_mean_vars)),
                        pearson_rho_low=rep(NA,length(exp1_mean_vars)), pearson_rho_high=rep(NA,length(exp1_mean_vars)), 
                        slope=rep(NA,length(exp1_mean_vars)), 
                        slope_low=rep(NA,length(exp1_mean_vars)), slope_high=rep(NA,length(exp1_mean_vars)))
  for (n in seq(1,length(exp1_mean_vars),1))  {
    # print(n)
    df_sub = subset(df_means, df_means$variable == exp1_mean_vars[n])
    
    
    if (length(unique((abs(df_sub$mean_value)))) > 1) {
      # Pearson correlation
      ci_pearson = ci_cor(data.frame(y1 = abs(df_sub[[indvar]])*-dir_to_agreement, 
                                     y2 = abs(df_sub$mean_value) ), 
                          method = "pearson", type = "normal", 
                          probs = c(alpha/length(exp1_mean_vars), 1-alpha/length(exp1_mean_vars)))
      df_mean_stat$pearson_rho[n] = ci_pearson[[3]]
      df_mean_stat$pearson_rho_low[n] = ci_pearson[[2]][1]
      df_mean_stat$pearson_rho_high[n] = ci_pearson[[2]][2]
    } else { 
    df_mean_stat$pearson_rho[n]    <- 0
    df_mean_stat$pearson_rho_low[n]  <- 0
    df_mean_stat$pearson_rho_high[n] <- 0
    }
    # LInear regression
    stat.lm <- lm(y2 ~ y1, data = tibble(y1 = abs(df_sub[[indvar]])*-dir_to_agreement,
                                         y2 = abs(df_sub$mean_value)))
    sd_slope <- summary(stat.lm)[[4]][2]
    df_mean_stat$slope[n] = stat.lm$coefficients[2]
    df_mean_stat$slope_low[n] = df_mean_stat$slope[n] - 1.96*sd_slope
    df_mean_stat$slope_high[n] = df_mean_stat$slope[n] + 1.96*sd_slope
  }
  df_mean_stat <- df_mean_stat[match(exp1_mean_vars, df_mean_stat$variable),]
  df_mean_stat$label <- factor(stats_labels, levels = stats_labels)
  df_mean_stat$is_pearson_rho_sig <-  !(0>df_mean_stat$pearson_rho_low & 0<df_mean_stat$pearson_rho_high)
  df_mean_stat$is_slope_sig <-  !(0>df_mean_stat$slope_low & 0<df_mean_stat$slope_high)
  
  
gg <- ggplot(data = df_mean_stat, aes(x = label, y = pearson_rho)) + 
  geom_point(shape=1, size=1) + ylim(-1,1) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_linerange(aes(ymin = pearson_rho_low, ymax = pearson_rho_high)) +
  ylab("Pearson's r") + xlab("Statistic") +
  scale_x_discrete(labels= parse(text = as.character(df_mean_stat$label))) +
  theme_classic(base_size=8) + theme(legend.position="none") 
print(gg)  
save_plot(paste(fig_path, fig_name, sep = ""), gg, ncol = 1, nrow = 1, 
          base_height = 1.75, base_asp = 3, base_width = 3, dpi = 600)


# Plot values of of mu, rmu, sigma, rsigma of d and b over simulations, and df
df_runs = tibble(Series = rep(seq(1,dim(df)[1],1),5),
                 param = c(rep("mu[DM]",dim(df)[1]), rep("r*mu[DM]",dim(df)[1]),
                 rep("sigma[pool]",dim(df)[1]), rep("r*sigma[pool]",dim(df)[1]),
                 rep("df[pool]",dim(df)[1])), 
                 value = c(df$mu_1dm, df$rmu_1dm, df$sigma_1d, df$rsigma_1d, df$df_1d)
                 )
df_runs$param <- factor(df_runs$param, levels = c("mu[DM]", "r*mu[DM]", "sigma[pool]","r*sigma[pool]","df[pool]"))


df_means <- df_runs %>% group_by(param) %>% summarize(Series=1,
  mean_value=mean(value),  is_constant = all(mean(value) == value),
  ymin = min(c(mean_value - 0.1 * mean_value, mean_value-.15)), #
  ymax = max(c(mean_value + 0.1 * mean_value, mean_value+.15))) #

# browser();


# df_means$ymin <-min(c(df_means$mean_value - 0.2 * df_means$mean_value, -.2))
# df_means$ymax <-max(c(df_means$mean_value + 0.2 * df_means$mean_value, +.2))

gg <- ggplot(data = df_runs, aes(x = Series, y = value)) +
  geom_line() +
  facet_wrap(vars(param), nrow=3,ncol=2,scales="free_y",labeller=label_parsed) +
  theme_classic(base_size=8) +
  theme(strip.text.x = element_text( margin = margin( b = 0, t = 0) )) 
  # geom_blank(data=df_means, aes(x = Series, y=mean_value, ymin = ymin, ymax = ymax))
gg
save_plot(paste(fig_path, str_replace(fig_name,"\\.[a-z]*$","_params.tiff"), sep = ""), gg, ncol = 1, nrow = 1, 
          base_height = 1.75, base_asp = 4, base_width = 3, dpi = 600)

save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
# load(file = "temp/debug.RData")


return(df_mean_stat)
}