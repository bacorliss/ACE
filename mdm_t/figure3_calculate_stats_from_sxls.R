

if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(readxl)
# Import table from Excel
source("R/row_stats_toolbox.R")



base_dir = "mdm_t"
fig_num = "3" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(getwd(),"/",base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE)


norm_points <- function(n) {x = rnorm(n,0,1); x=(x-mean(x))/sd(x); x}

parse_fract <- function(s) {
  #' Safe function to parse fraction from string instead of using eval
  x = strsplit(s,split="/")
  if (length(x[[1]])==2) {y=as.numeric(x[[1]][1])/as.numeric(x[[1]][2])
  } else {y=as.numeric(x)}; return(y)
}



calc_stats <- function(df, out_path, csv_name, p2_relrange) {
  #' Compute stats from data frame, one from each row, based on mean,std, and n 
  #' from each group
  #' #
  # browser();
  df_stat = data.frame(
                      x = factor(seq(1,dim(df)[1], 1)),
                      rmdm = rep(0,dim(df)[1]), 
                      mdm_ov_xbar = rep(0,dim(df)[1]),
                      mdm = rep(0,dim(df)[1]), 
                      cd = rep(0,dim(df)[1]),
                      p_nhst = rep(0,dim(df)[1]), 
                      p_equiv = rep(0,dim(df)[1]), 
                      p_2nd = rep(0,dim(df)[1]), 
                      bf = rep(0,dim(df)[1]),
                      xbar_dm = rep(0,dim(df)[1]), rxbar_dm = rep(0,dim(df)[1]),
                      wtt_lo = rep(0,dim(df)[1]),  wtt_hi = rep(0,dim(df)[1]), 
                      wtt_rlo = rep(0,dim(df)[1]), wtt_rhi = rep(0,dim(df)[1]))
  xbarA_str = paste("x","\U00AF","A",sep="")
  xbarB_str = paste("x","\U00AF","B",sep="")
  alphadm_str = paste("\U03B1","DM",sep="")
  for (n in seq(1, dim(df)[1]) ) {
    # Simulate data points for input into functions
    # Create normalized data with n1 points, times s1 plus x1
    xa = matrix(norm_points(as.double(df$nA[n])), nrow=1) * 
      as.double(df$sA[n]) + as.double(df[[xbarA_str]][n])
    # Create normalized data with n2 points, times s1 plus x2
    xb = matrix(norm_points(as.double(df$nB[n])), nrow=1) * 
      as.double(df$sB[n]) + as.double(df[[xbarB_str]][n])
    # Note: even though each xa,xb is one sample, we still caste them as matrices 
    # so we can use the row_effect_size functions to compute them (so it's the 
    # same function call as the rest of the paper)
    
    # Input points into each function
    df_stat$bf[n] <- row_bayesf_2s(xa, xb)
    df_stat$p_nhst[n] <- t.test(xa, xb, conf.level = 
                                  1 - parse_fract(df[[alphadm_str]][n]))$p.value
    df_stat$mdm[n] <- row_mdm(xa, xb, conf.level = 
                                        1 - parse_fract(df[[alphadm_str]][n]))
    df_stat$rmdm[n] <- 
      row_rmdm(m_c = xa, m_e = xb, conf.level = 1 - parse_fract(df[[alphadm_str]][n]))
  
    df_stat$mdm_ov_xbar[n] <- df_stat$mdm[n] / as.double(df[[xbarA_str]][n])
    
    df_stat$p_equiv[n] <- row_tost_2s(xa, xb, conf.level = 
                                        1 - parse_fract(df[[alphadm_str]][n]))
    df_stat$cd[n] <- row_cohend(xb, xa)
    # df_stat$p_2nd[n] <- row_sgpv(xa, xb, -1,1)
    df_stat$p_2nd[n] <- row_sgpv(xa, xb, -p2_relrange * as.double(df[[xbarA_str]][n]),
                                 p2_relrange * as.double(df[[xbarA_str]][n]))
    
    
    wtt <- t.test(xb,xa, var.equal = FALSE, conf.level = 1 - parse_fract(df[[alphadm_str]][n]))
    df_stat$wtt_lo[n]  <- wtt$conf.int[1]
    df_stat$wtt_hi[n]  <- wtt$conf.int[2]
    df_stat$xbar_dm[n] <- (mean(xb) - mean(xa))
    
    df_stat$wtt_rlo[n]  <- wtt$conf.int[1]/mean(xa)
    df_stat$wtt_rhi[n]  <- wtt$conf.int[2]/mean(xa)
    df_stat$rxbar_dm[n] <- (mean(xb) - mean(xa))/mean(xa)
  }
 
  
  # plot relative confidence intervals
  gg = ggplot(data = df_stat, aes(x=x, y = rxbar_dm)) +
    geom_pointrange(size = .4, aes(ymin = wtt_rlo, ymax = wtt_rhi)) + 
    geom_hline(yintercept=0, size = .5) +
   
    ylab( parse(text=paste("Relative~CI~of~mu[DM]"))) + xlab("Study #") +
    theme_classic(base_size=8)
  print(gg)
  save_plot(paste(out_path, "/rel_ci_", csv_name, ".tiff",sep=""),
            gg, ncol = 1, nrow = 1, base_height = 1.75, base_width = 6.5, dpi = 600)
  
  
  # Export results to disk
  write.csv(df_stat, paste(out_path,"/", csv_name,sep=""))
  
  return(df_stat)
}




# Total Cholesterol, Null Results
raw_chol_null <- read_excel("mdm_t/atherosclerosis_review.xlsx", sheet = "Chol Null", skip = 1)
chol_null <- raw_chol_null[-seq(min(which(is.na(raw_chol_null$sA))), dim(raw_chol_null)[1], 1),]
df_stat_chol_null <- calc_stats(chol_null, fig_path,"chol_null.csv", p2_relrange = 0.3)


# Total CHolesterol, crit Results
raw_chol_crit <- read_excel("mdm_t/atherosclerosis_review.xlsx", sheet = "Chol Crit", skip=1)
chol_crit <- raw_chol_crit[-seq(min(which(is.na(raw_chol_crit$sA))), dim(raw_chol_crit)[1], 1),]
df_stat_chol_crit <- calc_stats(chol_crit, fig_path,"chol_crit.csv", p2_relrange = 0.3)



# Plaque Size, Null Results
raw_plaq_null <- read_excel("mdm_t/atherosclerosis_review.xlsx", sheet = "Plaque Null", skip=1)
plaq_null <- raw_plaq_null[-seq(min(which(is.na(raw_plaq_null$sA))), dim(raw_plaq_null)[1], 1),]
df_stat_plaq_null <- calc_stats(plaq_null, fig_path,"plaq_null.csv", p2_relrange = 0.4)


# Plaque Size, Null Results
raw_plaq_crit <- read_excel("mdm_t/atherosclerosis_review.xlsx", sheet = "Plaque Crit", skip=1)
plaq_crit <- raw_plaq_crit[-seq(min(which(is.na(raw_plaq_crit$sA))), dim(raw_plaq_crit)[1], 1),]
df_stat_plaq_crit <- calc_stats(plaq_crit, fig_path,"plaq_crit.csv", p2_relrange = 0.4)
