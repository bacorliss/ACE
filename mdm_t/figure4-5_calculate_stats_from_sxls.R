

source("R/candidate_stats_from_xlsx.R")


base_dir = "mdm_t"
fig_num = "3" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(getwd(),"/",base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE)


# Total Cholesterol, Null Results
df_stat_chol_null <- compute_stats_from_xlsx("mdm_t/applied_examples.xlsx", "Chol Null", 
                                             fig_path,"chol_null.csv", rel_delta = 0.3)

# Total CHolesterol, Crit Results
df_stat_chol_crit <- compute_stats_from_xlsx("mdm_t/applied_examples.xlsx", "Chol Crit",
                                             fig_path,"chol_crit.csv", rel_delta = 0.3)

# Plaque Size, Null Results
df_stat_plaq_null <- compute_stats_from_xlsx("mdm_t/applied_examples.xlsx","Plaque Null", 
                                             fig_path,"plaq_null.csv", rel_delta = 0.4)

# Plaque Size, Crit Results
df_stat_plaq_crit <- compute_stats_from_xlsx("mdm_t/applied_examples.xlsx","Plaque Crit",
                                             fig_path,"plaq_crit.csv", rel_delta = 0.4)







