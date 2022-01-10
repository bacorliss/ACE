


# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
base_dir = "mdm_t"

#' Take same simulations fro similarity contests, but sweep the indepdent variable 
#' and plot error rate with 95% CI for subset of stats across indepdendent variable



# Calculate weighted power score for null results and positive results
#     For null results, weight each simulation by the PDF of at that mu/sigma, 
# scaled across the independent variable
#     For positive results, weight each simulation by 1-the PDF of at that mu/sigma, 
# scaled across indepdent variables