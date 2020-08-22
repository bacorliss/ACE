### Effect Size Contest



# TODO: for specifying offset distribution in generateData(), should the standard
# deviation for the control group be added by stds or variances?

library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(tidyr)
library(cowplot)
library(dplyr)
# Parallel processing
library(boot)
library(foreach)
library(doParallel)
library(stringr)
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
effect_size_dict[[2]] <- c("xdbar", "rxdbar", "sd", "rsd", "bf", "p_value",
                           "tostp", "cohen_d", "mmd",  "rmmd","nrand")
effect_size_dict[[3]] <- c("d2gtd1","2m1")
effect_size_dict[[4]] <- c("bar(x)", "r*bar(x)", "s", "r*s","Bf", "p[NHST]*phantom(.)",
                           "~p[TOST]", "Cd" , "delta[M]",
                           "r*delta[M]","Rnd")


pop_params_from_aoffset <- function( n_samples, n_obs, n_sims, 
                                    mus_1a, sigmas_1a, 
                                    mus_2a, sigmas_2a,
                                    mus_1ao, sigmas_1ao, 
                                    mus_2ao, sigmas_2ao) {
  # browser()
  # Calculate D based on A and offset from A parameters
  mus_1d = mus_1ao;  sigmas_1d = sigmas_1a + sigmas_1ao
  mus_2d = mus_2ao;  sigmas_2d = sigmas_2a + sigmas_2ao
  # calculate B from A and D
  mus_1b = mus_1a + mus_1d;  sigmas_1b = sqrt(sigmas_1d^2 - sigmas_1a^2)
  mus_2b = mus_2a + mus_2d;  sigmas_2b = sqrt(sigmas_2d^2 - sigmas_2a^2)
  # Initialize df with mu_1a, sigma_1a, mu_1b, sigma_1b, mu_1d,
  df = tibble( n_obs = n_obs, n_samples = n_samples,
    mu_1a = mus_1a, mu_1b = mus_1b, mu_1d = mus_1d,
    mu_2a = mus_2a, mu_2b = mus_2b, mu_2d = mus_2d, 
    sigma_1a = sigmas_1a, sigma_1b = sigmas_1b, sigma_1d = sigmas_1d,
    sigma_2a = sigmas_2a, sigma_2b = sigmas_2b, sigma_2d = sigmas_2d,
  )
  return(df)
}


pop_params_from_ab <- function( n_samples, n_obs, n_sims, 
                                mus_1a, sigmas_1a, 
                                mus_2a, sigmas_2a,
                                mus_1b, sigmas_1b, 
                                mus_2b, sigmas_2b) {
  # Calculate D based on A and B parameters
  mus_1d = mus_1b - mus_1a;  sigmas_1d = sqrt(sigmas_1a^2 + sigmas_1b^2)
  mus_2d = mus_2b - mus_2a;  sigmas_2d = sqrt(sigmas_2a^2 + sigmas_2b^2)
  # Initialize df with mu_1a, sigma_1a, mu_1b, sigma_1b, mu_1d,
  df = tibble( n_obs = n_obs, n_samples = n_samples,
               mu_1a = mus_1a, mu_2a = mus_2a, sigma_1a = sigmas_1a, sigma_2a = sigmas_2a,
               mu_1b = mus_1b, mu_2b = mus_2b, sigma_1b = sigmas_1b, sigma_2b = sigmas_2b,
               mu_1d = mus_1d, mu_2d = mus_2d, sigma_1d = sigmas_1d, sigma_2d = sigmas_2d,
  )
  return(df)
}


pop_params_switches <- function(df_init, switch_sign_mean_d, switch_sign_mean_ab, 
                                switch_group_ab, switch_exp_12) {
  
  
  df <- df_init
  # browser()
  
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
generateExperiment_Data <- function(n_samples, n_obs, n_sims, rand.seed,
                                    # Control group pop. parameters
                                    mus_1a, sigmas_1a, 
                                    mus_2a, sigmas_2a,
                                    # Experiment group pop. parameters
                                    mus_1b = NA, sigmas_1b = NA, 
                                    mus_2b = NA, sigmas_2b = NA,
                                    # Difference distribution pop. parameters
                                    mus_1ao = NA, sigmas_1ao = NA, 
                                    mus_2ao = NA, sigmas_2ao = NA,
                                    switch_group_ab = FALSE,
                                    switch_sign_mean_ab = FALSE,
                                    switch_sign_mean_d = FALSE,
                                    switch_exp_12 = FALSE,
                                    fig_name = "test.tiff",
                                    fig_path = "Figure/",
                                    label_dict = effect_size_dict,
                                    gt_colnames) {
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
  
  # Expand any singleton pop param arguments replicate to number of simulations
  input_args <- formalArgs(generateExperiment_Data)
  pargs <-grep("^(mus)|(sigmas)", input_args, value=TRUE)
  # For any pop param not equal in length to n_sims, expand
  for (n in seq_along(pargs)) {
    if (length(get(pargs[n]))==1) assign(pargs[n], rep(get(pargs[n]),n_sims))
  }

  # Record some parameter values for simulations
  set.seed(rand.seed +1)
  
  # Generate initial dataframe from params, no switching done yet  
  if (is.na(mus_1b)) {
    df_init <- pop_params_from_aoffset( n_samples, n_obs, n_sims, 
                                         mus_1a, sigmas_1a,  mus_2a, sigmas_2a,
                                         mus_1ao, sigmas_1ao, mus_2ao, sigmas_2ao) 
  } else {
    df_init   <- pop_params_from_ab( n_samples, n_obs, n_sims, 
                                   mus_1a, sigmas_1a,  mus_2a, sigmas_2a,
                                   mus_1b, sigmas_1b, mus_2b, sigmas_2b)
  }
  # Switch params if needed
  df <- pop_params_switches(df = df_init, switch_sign_mean_d = switch_sign_mean_d, 
                            switch_sign_mean_ab = switch_sign_mean_ab, 
                            switch_group_ab = switch_group_ab,
                            switch_exp_12 = switch_exp_12)
  
  # browser()
  
  # Define difference distribution    
  df$mu_1d <- df$mu_1b - df$mu_1a
  df$mu_2d <- df$mu_2b - df$mu_2a
  df$sigma_1d <- sqrt(df$sigma_1a^2 + df$sigma_1b^2)
  df$sigma_2d <- sqrt(df$sigma_2a^2 + df$sigma_2b^2) 
  
  # Calculate parameter of difference in means distribution (taken from mean of 
  # D since we had the option to invert the sign for D)
  df$mu_1md <- df$mu_1d
  df$mu_2md <- df$mu_2d
  df$sigma_1md <- df$sigma_1d /sqrt(n_obs)
  df$sigma_2md <- df$sigma_2d /sqrt(n_obs)
  
  # Calculate ratio of sigma_md/mu_md to determine how close D is close to zero,
  # and how absolute value folding will effect distribution.
  df$mu_ov_sigma_1md <- df$mu_1md / df$sigma_1md
  df$mu_ov_sigma_2md <- df$mu_2md / df$sigma_2md

  # Is: Exp2 mu[d] > Exp1 mu[d]
  df$is_mud_md2gtmd1 <-  abs(df$mu_2md) > abs(df$mu_1md)
  # Statistics of difference of means distribution 
  df$rmu_1md <- df$mu_1md/df$mu_1a
  df$rmu_2md <- df$mu_2md/df$mu_2a
  # Is: Exp2 relative_mu[d] > Exp1 relative_mu[d]
  df$is_rmud_md2gtmd1 <-  abs(df$rmu_2md) > abs(df$rmu_1md)
  
  # Is: pop_std2 > pop_std1 (TRUE)
  df$is_sigma_md2gtmd1 <-  df$sigma_2md > df$sigma_1md
  # Relative change variance of the mean compared to difference in means
  #   Realtive variance across scales
  # Is: rel_pop_std2 > rel_pop_std1, a/k/a CV1 > CV2
  df$rsigma_1md <- df$sigma_1md / abs(df$mu_1a + df$mu_1d/2)
  df$rsigma_2md <- df$sigma_2md / abs(df$mu_2a + df$mu_2d/2)
  df$is_rsigma_md2gtmd1 <-  df$rsigma_2md > df$rsigma_1md

  
  # Diff:  pop_mean2 - pop_mean1
  df$mean_mud_d2md1 <- df$mu_2md - df$mu_1md
  df$mean_rmud_d2md1 <- (df$mu_2md/df$mu_2a) - (df$mu_1md/df$mu_1a)
  # Diff:  pop_std2 - pop_std1
  df$mean_sigma_d2md1 <- df$sigma_2md - df$sigma_1md
  df$mean_rsigma_d2md1 <- df$sigma_2md/df$mu_2a - df$sigma_1md/df$mu_1a
  
  # Append columns for effect sizes, since multiple columns are used to analyze
  # each effect size, a dictionary of prefix, base, and suffix variable names 
  # are used.
  df[ paste(effect_size_dict$prefix[1], effect_size_dict$base,
            effect_size_dict$suffix[1], sep="_") ] <- rep(NaN,n_sims)
  df[ paste(effect_size_dict$prefix[2], effect_size_dict$base,
            effect_size_dict$suffix[2], sep="_") ] <- rep(NaN,n_sims)
  
  # Plot generated population parameters
  plot_population_params(df, fig_name = fig_name, fig_path = fig_path, 
                         gt_colnames = gt_colnames)
  
  
  # browser()
  return(df)
}


plot_population_params <- function(df_init, gt_colnames,fig_name,fig_path){
  
  # Output csv of agreement of input parameters to each individual input parameter
  param_fields = c("is_mud_md2gtmd1","is_rmud_md2gtmd1","is_sigma_md2gtmd1",
                   "is_rsigma_md2gtmd1")
  
  # Calculate indices of colname
  gt_param_inds <- unname(sapply(gt_colnames,function(x){pmatch(x,param_fields)}))
  gt_param_labels <- c("abs(~mu[DM]*phantom(.))", 
                       "abs(~r*mu[DM]*phantom(.))",
                       "~sigma[DM]*phantom(.)", "r*sigma[DM]")
  
  
  n_agreement = matrix(0, ncol = length(param_fields), nrow = length(param_fields))
  pwise_binom_p <- n_agreement
  params <- str_match(param_fields, "is_(.*)_md2gtmd1")[,2]
  
  for (r in seq(1,length(param_fields),1)) {
    for (c in seq(1,length(param_fields),1)) {
      n_agreement[r,c]   <- sum(df_init[[param_fields[r]]] == df_init[[param_fields[c]]])
      pwise_binom_p[r,c] <-
        prop.test(n_agreement[r,c], dim(df_init)[1], alternative = "two.sided",
                  conf.level = 1-0.05/(4*length(gt_colnames)), correct = TRUE)$p.value
    }
  }
  # Corrected p values based on number of independent variables selected * 4
  pwise_binom_p_corr <- pmin(4*length(gt_colnames) * pwise_binom_p,rep(1,prod(dim(pwise_binom_p))))
  str_binom_p <- matrix(sapply(pwise_binom_p_corr, function(x) if(x> 0.05) 
  {sprintf("%0.2f",x)} else {sprintf("%0.2e",x)}),nrow = dim(pwise_binom_p)[1])
  colnames(str_binom_p) <- params
  rownames(str_binom_p) <- params
  # Write to csv
  csv_path <- paste(fig_path, "params_",str_replace(fig_name, ".tiff$", ".csv"), sep="")
  cat("Table 1: Binomial test of agreement by group\n", file = csv_path)
  suppressWarnings(write.table(str_binom_p, csv_path, append = FALSE, 
                               col.names = TRUE, sep=","))
  
  # Calculate binomial confidence intervals for each sucess rate for each parameter
  binom_mu <- prop.test(sum(df_init$is_mud_md2gtmd1), dim(df_init)[1], 
                        conf.level=1-0.05/(4*length(gt_colnames)), correct = FALSE)
  binom_rmu <- prop.test(sum(df_init$is_rmud_md2gtmd1), dim(df_init)[1], 
                         conf.level=1-0.05/(4*length(gt_colnames)), correct = FALSE)
  binom_sigma <- prop.test(sum(df_init$is_sigma_md2gtmd1), dim(df_init)[1], 
                           conf.level=1-0.05/(4*length(gt_colnames)), correct = FALSE)
  binom_rsigma <- prop.test(sum(df_init$is_rsigma_md2gtmd1), dim(df_init)[1], 
                            conf.level=1-0.05/(4*length(gt_colnames)), correct = FALSE)
  df_params <- rbind(
    tibble(group="mu", estimate = binom_mu$estimate, 
           lci = binom_mu$conf.int[1], uci = binom_mu$conf.int[2]),
    tibble(group="rmu", estimate = binom_rmu$estimate, 
           lci = binom_rmu$conf.int[1], uci = binom_rmu$conf.int[2]),
    tibble(group="sigma", estimate = binom_sigma$estimate, 
           lci = binom_sigma$conf.int[1], uci = binom_sigma$conf.int[2]),
    tibble(group="rsigma", estimate = binom_rsigma$estimate, 
           lci = binom_rsigma$conf.int[1], uci = binom_rsigma$conf.int[2]))
  df_params$group <- factor(df_params$group, levels = df_params$group)
  
  # Plot confidence interval of success rate for each parameter and their agreement
  p <- ggplot(df_params, aes(x = group,  y = estimate)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey",linetype="dashed") +
    geom_linerange(aes(ymin = lci, ymax = uci), size = 0.5) +
    geom_point(size = 1.25, fill = "white", shape = 1) + 
    ylab("Fract Exp 1 < Exp 2    ") +
    xlab("Groundtruth") +
    theme_classic(base_size = 8) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(y=c(0,1), clip = "off") +
    scale_x_discrete(expand=c(0.1, 0), labels =
       c('mu' = parse(text=paste(gt_param_labels[1],"*phantom(.)")),
       'rmu'   = parse(text=paste("phantom(.)*",gt_param_labels[2])),
       'sigma' = parse(text=gt_param_labels[3]),
       'rsigma'   = parse(text=gt_param_labels[4])))
  # Add statistical annotation above plot, one for each ground truth variable 
  # specified (up to 2)
  p <- p +  theme(plot.margin = unit(c(13,3,3,3), "pt")) +
    annotate("text",label = paste(gt_param_labels[gt_param_inds[1]],"~phantom(.)",
                                  sep = ""), x = 0.8, size=2.5,
             y = 1.13,vjust = 0, hjust=1,parse=TRUE) +
    geom_text(y = 1.13, aes(label = ifelse(pwise_binom_p_corr[gt_param_inds[1],]<0.05, "#","")),
              size = 2.5, vjust=0, hjust=0.5) +
    annotate("segment", x = 0.8, xend = 4, y = 1.09, yend = 1.09, colour = "black", size=.2) 
  if (length(gt_colnames)==2){
    p <- p +  theme(plot.margin = unit(c(22,0,0,0), "pt")) +
      annotate("text",label = paste(gt_param_labels[gt_param_inds[2]],"~phantom(.)",sep=""),
               x = 0.8, size=2.5, y = 1.3,vjust = 0, hjust=1,parse=TRUE) +
      geom_text(y = 1.3,aes(label = ifelse(pwise_binom_p_corr[gt_param_inds[2],]<0.05, "#","")),
                size = 2.5, vjust=0, hjust=0.5) +
      annotate("segment", x = 0.8, xend = 4, y = 1.27, yend = 1.27, colour = "black", size=.2) 
  }    
  print(p)
  save_plot(paste(fig_path, 'gt_',fig_name, sep = ""), p, ncol = 1, nrow = 1, 
            base_height = 1.5, base_asp = 3, base_width = 1.35, dpi = 600)
  
  
  
  # Export csv file for agreement between each variable to others
  # Plot histogram of mu[D]/sigma[D] to demonstrate how far from zero D is  
  df <-tibble(group = as.factor(c(rep(1,dim(df_init)[1]),rep(2,dim(df_init)[1]))),
              mu_ov_sigma = c(df_init$mu_1md/df_init$sigma_1md,
                              df_init$mu_2md/df_init$sigma_2md))
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
  # browser();
}



quantify_esize_simulation <- function(df, include_bf = FALSE, rand.seed = 0, 
                                      parallelize_bf = FALSE) {
  if (dim(df)[1] != 1) stop("Need to specify a single row of effect size matrix")
  
  set.seed(rand.seed)
  # Transform simulated samples with normal parameters (mean and std) for 
  # current round of simulation
  
  # Use Exp 1 and 2 coefficients to generate data from normalized base data
  x_1a = matrix(rnorm(df$n_samples * df$n_obs, mean = df$mu_1a, 
                      sd = df$sigma_1a), nrow = df$n_samples, 
                ncol = df$n_obs)
  x_1b = matrix(rnorm(df$n_samples * df$n_obs, mean = df$mu_1b, 
                      sd = df$sigma_1b), nrow = df$n_samples, 
                ncol = df$n_obs)
  
  x_2a = matrix(rnorm(df$n_samples * df$n_obs, mean = df$mu_2a, 
                      sd = df$sigma_2a), nrow = df$n_samples, 
                ncol = df$n_obs)
  x_2b = matrix(rnorm(df$n_samples * df$n_obs, mean = df$mu_2b, 
                      sd = df$sigma_2b), nrow = df$n_samples, 
                ncol = df$n_obs)
  
  # Sample estimate of difference in means
  xbar_1d = rowMeans(x_1b) - rowMeans(x_1a)
  xbar_2d = rowMeans(x_2b) - rowMeans(x_2a)
  
  # Sample estimate of standard deviation of difference in means
  s_1md = sqrt(rowSds(x_1a)^2/df$n_obs + rowSds(x_1b)^2/df$n_obs)
  s_2md = sqrt(rowSds(x_2a)^2/df$n_obs + rowSds(x_2b)^2/df$n_obs)
  
  # Basic Summary statistical comparisons
  # Means
  diff_xdbar = abs(xbar_2d) - abs(xbar_1d)
  df$fract_xdbar_d2gtd1 = sum(diff_xdbar > 0) / df$n_samples
  df$mean_diff_xdbar_2m1  = mean(diff_xdbar)
  
  # Stds
  diff_sd = s_2md - s_1md
  df$fract_sd_d2gtd1  = sum(diff_sd > 0) / df$n_samples
  df$mean_diff_sd_2m1 = mean(diff_sd)
  
  # Rel Means: mean divided by control mean
  diff_rxdbar = abs(xbar_2d/rowMeans(x_2a)) - abs(xbar_1d/rowMeans(x_1a))
  df$fract_rxdbar_d2gtd1 = sum( diff_rxdbar > 0) / df$n_samples
  df$mean_diff_rxdbar_2m1  =  mean(diff_rxdbar)
  
  # Rel STDs: sd divided by control mean
  diff_rsd = s_2md / (rowMeans(x_2a) + xbar_2d/2) -
    s_1md / (rowMeans(x_1a) + xbar_1d/2)
  df$fract_rsd_d2gtd1 = sum( diff_rsd > 0) / df$n_samples
  df$mean_diff_rsd_2m1  = mean(diff_rsd)
  
  if (include_bf) {
    # Bayes factor
    bf_1 = row_bayesf_2s(x_1a, x_1b, parallelize = parallelize_bf, paired = FALSE)
    bf_2 = row_bayesf_2s(x_2a, x_2b, parallelize = parallelize_bf, paired = FALSE)
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
  p_value_1 = 2*pnorm(-abs(z_score_1))
  p_value_2 = 2*pnorm(-abs(z_score_2))
  diff_p_value <-  p_value_2 - p_value_1
  # Must switch signs
  df$fract_p_value_d2gtd1 = sum(-diff_p_value > 0) / df$n_samples
  df$mean_diff_p_value_2m1 = mean(diff_p_value)
  
  # TOST p value (Two tailed equivalence test)
  tostp_1 <- row_tost_2s(x_1b, x_1a,low_eqbound = -.1,high_eqbound = .1)
  tostp_2 <- row_tost_2s(x_2b, x_2a,low_eqbound = -.1,high_eqbound = .1)
  diff_tostp <- tostp_2 - tostp_1
  df$fract_tostp_d2gtd1 <- sum(diff_tostp > 0) / df$n_samples
  df$mean_diff_tostp_2m1 <- mean(diff_tostp) 
  
  # Delta Family of effect size
  # Cohens D
  diff_cohen_d = abs(row_cohend(x_2a, x_2b)) - abs(row_cohend(x_1a, x_1b))
  df$fract_cohen_d_d2gtd1 = sum(diff_cohen_d > 0) / df$n_samples
  df$mean_diff_cohen_d_2m1 =  mean(diff_cohen_d)
  
  # Most Mean Diff
  mmd_1 = row_mmd_2s(x_1a, x_1b, paired = FALSE)
  mmd_2 = row_mmd_2s(x_2a, x_2b, paired = FALSE)
  diff_most_mean_diff = mmd_2 - mmd_1
  df$fract_mmd_d2gtd1 = sum(diff_most_mean_diff > 0) / df$n_samples
  df$mean_diff_mmd_2m1 = mean(diff_most_mean_diff)
  
  # Relative Most Mean Diff
  diff_rmmd = mmd_2 / rowMeans(x_2a) -
    mmd_1 / rowMeans(x_1a)
  df$fract_rmmd_d2gtd1 = sum(diff_rmmd > 0) / df$n_samples
  df$mean_diff_rmmd_2m1 = mean(diff_rmmd)
  
  
  # Random group
  diff_nrand = rowMeans(matrix(rnorm(df$n_samples * df$n_obs, mean = 0, sd = 1), 
                               nrow = df$n_samples, ncol = df$n_obs))
  df$fract_nrand_d2gtd1 = sum(diff_nrand > 0 ) / df$n_samples
  df$mean_diff_nrand_2m1 = mean(diff_nrand)
    
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
  
  n_sims = dim(df_in)[1]
  # Initialize output data frame
  df <- df_in
  
  # Only perform simulations if results not saved to disk
  if (!file.exists(paste(out_path,data_file_name,sep="")) | overwrite) {
    
    # browser();
    if (parallel_sims) {
      print("starting parallel processing")
      #setup parallel back end to use many processors
      cores = detectCores()
      cl = makeCluster(cores[1]-1)
      registerDoParallel(cl)
      
      df <- foreach(n = 1:n_sims, .combine = rbind, .packages = 
                      c("BayesFactor","TOSTER"),
                .export = c(row_effect_sizes_fun, mmd_functions,
                            "quantify_esize_simulation")) %dopar% {
        #calling a function
        tempMatrix = quantify_esize_simulation(df[n,], include_bf, 
                                               rand.seed = rand.seed+n,
                                               parallelize_bf = FALSE) 
        tempMatrix
      }
      #stop cluster
      stopCluster(cl)
    }else{
      # Process effect sizes serially
      for (n in seq(1,n_sims,1)) {
        df[n,] <- quantify_esize_simulation(df[n,], include_bf, rand.seed = rand.seed+n, 
                                            parallelize_bf = FALSE) 
      }
    }
    
    # Save dataframed results to a file
    saveRDS(df, file = paste(out_path,data_file_name,sep=""))
  } else {
    # Restore the data frame results from disk
    df <- readRDS(file = paste(out_path,data_file_name,sep=""))
  }
  # browser()
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


plot_esize_simulations <- function(df_pretty, fig_name, fig_path, y_ax_str) {
  
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

  # Basic violin plot
  p <- ggplot(df_result, aes(x=name,  y=mean, group=name)) +
    geom_hline(yintercept = 0.5, size=0.5, color="grey") +
    geom_linerange(aes(ymin = bs_ci_mean_lower, ymax = bs_ci_mean_upper), size = 0.5) +
    geom_point(size=1,fill="white", shape = 1) + 
    xlab("Statistic") +
    ylab(parse(text=paste("Error~Rate~Lower~phantom(.)*", y_ax_str,"~phantom(.)"))) +
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

process_esize_simulations <- function(df_init, gt_colname, y_ax_str, out_path="temp/",
                                      fig_name,fig_path,var_suffix = "fract",include_bf = TRUE,
                                      parallel_sims = TRUE) {
  
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
  df_plotted <- plot_esize_simulations(df = df_pretty, fig_name = fig_name, 
                                       fig_path = fig_path, y_ax_str = y_ax_str)
  
  
  all_dfs <- vector(mode="list", length=4)
  names(all_dfs) <- c("df_es", "df_tidy", "df_pretty", "df_plotted")
  all_dfs[[1]] <- df_es; all_dfs[[2]] <- df_tidy; 
  all_dfs[[3]] <- df_pretty; all_dfs[[4]] <- df_plotted; 
  return(all_dfs)
}