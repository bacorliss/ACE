
#' Produce stats look-up table for calculating MDM


# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(broom)
p_load(scales)
p_load(ggplot2)
p_load(dplyr)
p_load(plyr)
p_load(grid)
p_load(gridExtra)
p_load(colorspace)
p_load("RColorBrewer")
p_load(cowplot)
# User defined libraries
source("R/mdm.R")
source("R/row_stats_toolbox.R")
base_dir = "mdm_t"

# Figure parameters
#-------------------------------------------------------------------------------
out_path = file.path(getwd(), paste(base_dir, "/figure/T1",sep=""))
dir.create(out_path, showWarnings = FALSE)

# Simulation parameters
#-------------------------------------------------------------------------------
mu_ov_sigma = c(seq(0, 0.5, 0.1), 0.75, 1, 2,3)
n_samples = 35
n_obs = 50
set.seed(0)


# Calculate probability by zscore
zf = seq(0.1, 4.5, 0.1)
fn = matrix(rep(0,length(mu_ov_sigma)*length(zf)), ncol = length(mu_ov_sigma), nrow = length(zf))
colnames(fn) <- mu_ov_sigma
rownames(fn) <- zf
for (r in seq_along(zf)) {
  for (c in seq_along(mu_ov_sigma)) {  
    fn[r,c] <- pnorm(zf[r], mu_ov_sigma[c], 1) - pnorm(-zf[r], mu_ov_sigma[c], 1)
  }
}
write.table(fn, file = file.path(out_path, "zf_area_folded_normal.csv"), append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, fileEncoding = "")



# Calculate CDF at fixed offsets from mu/sigma
dzf = c(seq(0, 0.45,0.05), seq(0.50, 4.1, 0.1))
# dzf = c(seq(-4.5, 4.5,0.1))

dfn = matrix(rep(0,length(mu_ov_sigma)*length(dzf)), ncol = length(mu_ov_sigma), nrow = length(dzf))
colnames(dfn) <- mu_ov_sigma
rownames(dfn) <- dzf

for (r in seq_along(dzf)) {
  for (c in seq_along(mu_ov_sigma)) {  
    dfn[r,c] <- pnorm(dzf[r]+mu_ov_sigma[c], mu_ov_sigma[c], 1) - 
      pnorm(-dzf[r]-mu_ov_sigma[c], mu_ov_sigma[c], 1)
  }
}
write.table(dfn, file = file.path(out_path, "dzf_area_folded_normal.csv"),
            append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, fileEncoding = "")

# Calculate CDF at fixed offsets from mu/sigma
cumprob = c(seq(0.5, 0.6, 0.02), seq(0.61, 0.99,0.01), .995, .999, .9995, .9999)
mu_ov_sigma = c(seq(0, 0.5, 0.1), 1, 2,3,6)
# dzf = c(seq(-4.5, 4.5,0.1))

zfn = matrix(rep(0,length(mu_ov_sigma)*length(cumprob)), ncol = length(mu_ov_sigma), 
             nrow = length(cumprob))
colnames(zfn) <- mu_ov_sigma
rownames(zfn) <- cumprob
fun_pnstandnorm <- function (x, mu, sigma) {pnorm(x, mean = mu, sd = sigma) -
    pnorm(-x, mu, sigma)}

for (r in seq_along(cumprob)) {
  for (c in seq_along(mu_ov_sigma)) {  
    zfn[r,c] <- 
      uniroot(function(x) fun_pnstandnorm(x, mu = mu_ov_sigma[c], sigma = 1) - cumprob[r],
      interval = c(0, mu_ov_sigma[c]+ 10*qnorm(p = max(cumprob))),  tol = .Machine$double.eps)$root - mu_ov_sigma[c]
      # uniroot(function(x) 
      #   pnorm( x + mu_ov_sigma[c], mean = mu_ov_sigma[c], sd = 1) - 
      #   pnorm(-x + mu_ov_sigma[c], mean = mu_ov_sigma[c], sd = 1) - cumprob[r], 
      #   interval = c(0,  10*qnorm(p=max(cumprob))),
      #   tol = .Machine$double.eps)$root
  }
}
write.table(zfn, file = file.path(out_path, "dcumprob_area_folded_normal.csv"),
            append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, fileEncoding = "")




