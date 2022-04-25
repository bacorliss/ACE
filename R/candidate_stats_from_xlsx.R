

if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(readxl)
p_load(ggplot2)
p_load(cowplot)
p_load(stringr)

# Import table from Excel
source("R/row_stats_toolbox.R")

norm_points <- function(n) {x = rnorm(n,0,1); x=(x-mean(x))/sd(x); x}

parse_fract <- function(s) {
  #' Safe function to parse fraction from string instead of using eval
  x = strsplit(s,split="/")
  if (length(x[[1]])==2) {y=as.numeric(x[[1]][1])/as.numeric(x[[1]][2])
  } else {y=as.numeric(x)}; return(y)
}

compute_stats_from_xlsx <- function(xlsx_path, xlsx_sheet_name, out_path, 
                                    out_csv_name, rel_delta,
                                    rel_delta_sign = c(-1,1), rylim=c(-1.5, 1.5), 
                                    rand.seed = 0, delta_color = "#00B050") {
  #' Compute stats from data frame, one from each row, based on mean,std, and n 
  #' from each group
  #' #
  # rel_delta_sign 3 way flag (0 both signs, -1, +1)
  # Load xlsx at specified sheet and remove empty rows
  raw_df <- read_excel(xlsx_path, sheet = xlsx_sheet_name, skip = 1)
  df <- raw_df[-seq(min(which(is.na(raw_df$sX))), dim(raw_df)[1], 1),]
  
  
  
  save(list = ls(all.names = TRUE), file = "temp/compute_stats_from_xlsx.RData",
       envir = environment())
  # load(file = "temp/compute_stats_from_xlsx.RData")
  
  
  df_stat = data.frame(
    x = factor(seq(1,dim(df)[1], 1)),
    rmdm = rep(0,dim(df)[1]), mdm = rep(0,dim(df)[1]), 
    rldm = rep(0,dim(df)[1]), ldm = rep(0,dim(df)[1]), 
    cd = rep(0,dim(df)[1]),
    p_nhst = rep(0,dim(df)[1]), 
    p_equiv = rep(0,dim(df)[1]), 
    p_2nd = rep(0,dim(df)[1]), 
    bf = rep(0,dim(df)[1]),
    xbar_dm = rep(0,dim(df)[1]), rxbar_dm = rep(0,dim(df)[1]),
    cred_bound_lo = rep(0,dim(df)[1]),  cred_bound_hi = rep(0,dim(df)[1]),
    cred_bound_lo_ovc = rep(0,dim(df)[1]),  cred_bound_hi_ovc = rep(0,dim(df)[1]), 
    rcred_bound_lo = rep(0,dim(df)[1]), rcred_bound_hi = rep(0,dim(df)[1])
    )
  
  xbar_str = paste("x","\U00AF",sep="")
  ybar_str = paste("y","\U00AF",sep="")
  alphadm_str = paste("\U03B1","DM",sep="")
  for (n in seq(1, dim(df)[1]) ) {
    # Simulate data points for input into functions
    # Create normalized data with n1 points, times s1 plus x1
    # browser()
    xa = matrix(norm_points(as.double(df$m[n])), nrow=1) * 
      as.double(df$sX[n]) + as.double(df[[xbar_str]][n])
    # Create normalized data with n2 points, times s1 plus x2
    xb = matrix(norm_points(as.double(df$n[n])), nrow=1) * 
      as.double(df$sY[n]) + as.double(df[[ybar_str]][n])
    # Note: even though each xa,xb is one sample, we still caste them as matrices 
    # so we can use the row_effect_size functions to compute them (so it's the 
    # same function call as the rest of the paper)
    
    pooled_var = sqrt(sd(xa)^2/length(xa) + sd(xb)^2/length(xb))
    alpha_dm = parse_fract(df[[alphadm_str]][n])
    
    # Input points into each function
    df_stat$bf[n] <- row_bayesf_2s(xa, xb, deltas = rel_delta*mean(xa)/pooled_var)
    df_stat$p_nhst[n] <- t.test(xa, xb, conf.level = 
                                  1 - alpha_dm)$p.value
    df_stat$mdm[n] <- row_mdm(xa, xb, conf.level = 
                                1 - alpha_dm, rand.seed = rand.seed)
    df_stat$rmdm[n] <- row_rmdm(m_c = xa, m_e = xb, conf.level = 
                                  1 - alpha_dm, rand.seed = rand.seed)
    
    df_stat$ldm[n] <- row_ldm(xa, xb, conf.level = 
                                1 - alpha_dm, rand.seed = rand.seed)
    df_stat$rldm[n] <- 
      row_rldm(m_c = xa, m_e = xb, conf.level = 
                 1 - alpha_dm, rand.seed = rand.seed)
    
    
    # df_stat$mdm_ov_xbar[n] <- df_stat$mdm[n] / as.double(df[[xbar_str]][n])
    
    df_stat$p_equiv[n] <- row_tost_2s(xa, xb, deltas = rel_delta*mean(xa), 
                                      conf.level = 1 - alpha_dm)
    df_stat$cd[n] <- row_cohend(xb, xa)
    
    df_stat$p_2nd[n] <- row_sgpv(xa, xb, -rel_delta * mean(xa), rel_delta * mean(xa),conf.level = 1 - alpha_dm)
    
    
    bounds <- credint(c(xa), c(xb), conf.level = 1 - alpha_dm, relative = FALSE, 
                      rand.seed = rand.seed)
    rbounds <- credint(c(xa), c(xb), conf.level = 1 - 2*alpha_dm, relative = TRUE, 
                       rand.seed = rand.seed)
    wtt <- t.test(xb,xa, var.equal = FALSE, conf.level = 1 - alpha_dm)
    df_stat$cred_bound_lo[n]  <- bounds[1]
    df_stat$cred_bound_hi[n]  <- bounds[2]
    df_stat$xbar_dm[n] <- (mean(xb) - mean(xa))
    
    
    df_stat$cred_bound_lo_ovc[n]  <- bounds[1]/mean(xa)
    df_stat$cred_bound_hi_ovc[n]  <- bounds[2]/mean(xa)
    
    df_stat$rcred_bound_lo[n]  <- rbounds[1]
    df_stat$rcred_bound_hi[n]  <- rbounds[2]
    df_stat$rxbar_dm[n] <- (mean(xb) - mean(xa))/mean(xa)
  }
  print(df_stat)
  
  # plot unscaled credible intervals
  gg = ggplot(data = df_stat, aes(x=x, y = xbar_dm)) +
    geom_pointrange(size = .75, fatten = 1.5, aes(ymin = cred_bound_lo, ymax = cred_bound_hi)) + 
    geom_hline(yintercept=0, size = .5) +
    ylab( parse(text=paste("Cred.~Int.~of~mu[DM]"))) + xlab("Study #") +
    theme_classic(base_size=8)
  print(gg)
  save_plot(paste(out_path, "/credint_", out_csv_name, ".tiff",sep=""),
            gg, ncol = 1, nrow = 1, base_height = 1.75, base_width = 6.5, dpi = 600)
  
  
  # plot relative confidence intervals
  gg = ggplot(data = df_stat, aes(x=x, y = rxbar_dm)) +
    geom_pointrange(size = .75, fatten = 1.5,  aes(ymin = rcred_bound_lo, ymax = rcred_bound_hi)) + 
    geom_hline(yintercept=0, size = .5) +
    geom_hline(yintercept = rel_delta * rel_delta_sign, color = delta_color) + 
    ylab( parse(text=paste("Cred.~Int.~of~r*mu[DM]"))) + xlab("Study #") +
    theme_classic(base_size=8) +
    coord_cartesian(ylim=rylim)
  print(gg)
  save_plot(paste(out_path, "/rcredint_", out_csv_name, ".tiff",sep=""),
            gg, ncol = 1, nrow = 1, base_height = 2, base_width = 6.5, dpi = 600)
  
  
  # Export results to disk
  write.csv(df_stat, paste(out_path,"/", out_csv_name,sep=""))
  
  write.csv(t(df_stat), paste(out_path,"/", str_replace(out_csv_name,'.csv','_t.csv'),sep=""))
  # browser();
  
  return(df_stat)
}

