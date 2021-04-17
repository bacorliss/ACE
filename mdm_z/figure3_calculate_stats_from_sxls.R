

if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(readxl)
# Import table from Excel
source("R/row_stats_toolbox.R")



base_dir = "mdm_z"
fig_num = "3" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(getwd(),"/",base_dir,"/figure/F",fig_num, sep="")
dir.create(file.path(getwd(), fig_path), showWarnings = FALSE)


norm_points <- function(n) {x = rnorm(n,0,1); x=(x-mean(x))/sd(x); x}

parse_fract <- function(s) {
  #' Safe function to parse fraction from string instead of using eval
  x = strsplit(s,split="/")
  if (length(x[[1]])==2) {y=as.numeric(x[[1]][1])/as.numeric(x[[1]][2])
  } else {y=as.numeric(x)}; return(y)
}



calc_stats <- function(df, out_path, csv_name) {
  #' Compute stats from data frame, one from each row, based on mean,std, and n 
  #' from each group
  #' #
  df_stat = data.frame(dm = rep(0,dim(df)[1]), rdm = rep(0,dim(df)[1]), cd = rep(0,dim(df)[1]),
                       p_nhst = rep(0,dim(df)[1]), p_equiv = rep(0,dim(df)[1]), 
                       bf = rep(0,dim(df)[1]))
  
  for (n in seq(1, dim(df)[1]) ) {
    # Simulate data points for input into functions
    # Create normalized data with n1 points, times s1 plus x1
    xa = matrix(norm_points(as.double(df$n1[n])), nrow=1) * 
      as.double(df$s1[n]) + as.double(df[["x̅1"]][n])
    # Create normalized data with n2 points, times s1 plus x2
    xb = matrix(norm_points(as.double(df$n2[n])), nrow=1) * 
      as.double(df$s2[n]) + as.double(df[["x̅2"]][n])
    
    # Input points into each function
    df_stat$bf[n] <- row_bayesf_2s(xa, xb)
    df_stat$p_nhst[n] <- t.test(xa, xb, conf.level = 1 - parse_fract(df[["αDM"]][n]))$p.value
    df_stat$dm[n] <- row_mdm_2s_zdist(xa, xb, conf.level = 1 - parse_fract(df[["αDM"]][n]))
    df_stat$rdm[n] <- df_stat$dm[n] / as.double(df[["x̅1"]][n])
    df_stat$p_equiv[n] <- row_tost_2s(xa, xb, conf.level = 1 - parse_fract(df[["αDM"]][n]))
    df_stat$cd[n] <- row_cohend(xa, xb)
  }
 
  # Export results to disk
  write.csv(df_stat, paste(out_path,"/", csv_name,sep=""))
  
  return(df_stat)
}




# Total CHolesterol, Null Results
raw_chol_null <- read_excel("mdm_z/Atherosclerosis_Review.xlsx", sheet = "Chol Null", skip=1)
chol_null <- raw_chol_null[-seq(min(which(is.na(raw_chol_null$s1))), dim(raw_chol_null)[1], 1),]
df_stat_chol_null <- calc_stats(chol_null, fig_path,"chol_null.csv")


# Total CHolesterol, crit Results
raw_chol_crit <- read_excel("mdm_z/Atherosclerosis_Review.xlsx", sheet = "Chol Crit", skip=1)
chol_crit <- raw_chol_crit[-seq(min(which(is.na(raw_chol_crit$s1))), dim(raw_chol_crit)[1], 1),]
df_stat_chol_crit <- calc_stats(chol_crit, fig_path,"chol_crit.csv")



# Plaque Size, Null Results
raw_plaq_null <- read_excel("mdm_z/Atherosclerosis_Review.xlsx", sheet = "Plaque Null", skip=1)
plaq_null <- raw_plaq_null[-seq(min(which(is.na(raw_plaq_null$s1))), dim(raw_plaq_null)[1], 1),]
df_stat_plaq_null <- calc_stats(plaq_null, fig_path,"plaq_null.csv")


# Plaque Size, Null Results
raw_plaq_crit <- read_excel("mdm_z/Atherosclerosis_Review.xlsx", sheet = "Plaque Crit", skip=1)
plaq_crit <- raw_plaq_crit[-seq(min(which(is.na(raw_plaq_crit$s1))), dim(raw_plaq_crit)[1], 1),]
df_stat_plaq_crit <- calc_stats(plaq_crit, fig_path,"plaq_crit.csv")

