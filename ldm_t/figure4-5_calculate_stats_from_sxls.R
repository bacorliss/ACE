


source("R/candidate_stats_from_xlsx.R")



base_dir = "ldm_t"
fig_num = "3" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(getwd(),"/",base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE)


# Total CHolesterol, crit Results
df_stat_chol_crit <- compute_stats_from_xlsx("ldm_t/applied_examples.xlsx", "Chol Crit", 
                                fig_path,"chol_crit.csv", rel_delta = 0.20, rel_delta_sign = -1, rylim = c(-1., 0),
                                delta_color = "#0070C0")



# Plaque Size, Null Results
df_stat_plaq_crit <- compute_stats_from_xlsx("ldm_t/applied_examples.xlsx","Plaque Crit", 
                                fig_path,"plaq_crit.csv", rel_delta = 0.25, rel_delta_sign = -1, rylim = c(-1.25,0),
                                delta_color = "#0070C0")


