


library(broom)
library(scales)
library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(colorspace)
library("RColorBrewer")
library(cowplot)
source("R/mmd.R")
source("R/row_effect_sizes.R")



dir.create(file.path(getwd(), paste("figure/T1",sep="")), showWarnings = FALSE)

mu_ov_sigma = c(seq(0, 0.5, 0.1), 0.75, 1, 2,3)



n_samples = 35
n_obs = 50
set.seed(0)

zf = seq(0.1, 4.5, 0.1)
fn = matrix(rep(0,length(mu_ov_sigma)*length(zf)), ncol = length(mu_ov_sigma), nrow = length(zf))
colnames(fn) <- mu_ov_sigma
rownames(fn) <- zf



for (r in seq_along(zf)) {
  for (c in seq_along(mu_ov_sigma)) {  
    fn[r,c] <- pnorm(zf[r], mu_ov_sigma[c], 1) - pnorm(-zf[r], mu_ov_sigma[c], 1)
  }
}
write.table(fn, file = "figure/T1/zf_area_folded_normal.csv", append = FALSE, quote = TRUE, sep = ",",
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
write.table(dfn, file = "figure/T1/dzf_area_folded_normal.csv", append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, fileEncoding = "")


# Export LU table to disk
# df_lut = data.frame(abs_nmu = df_coeff$mu, coeff_mmd_95 = df_coeff$coeff_mmd_95)
# write.csv(x=df_lut, file=file.path(getwd(),"figure/T2/coeff_mmd_CLa_CL2a.csv"))
