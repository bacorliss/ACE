

# Load packages
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)
library(rlang)
library(colorspace)
library(VGAM)
library(boot)
library(dplyr)
library(cowplot)



library(VGAM)


fig_num = "2"
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)

# dnorm: density function: proability for a specific value of x/zscore of PDF
# pnorm: cum. probability/quantile for a given zscore
# qnorm: return z_score for a given cum. probability/quantile



tri_color <- brewer.pal(n = 3, name = "Set1")
distrs=c("FN", "NN","SN")

# Subfigure A: Compare FN,NN,SN with mu = 1 and swept sigma
#------------------------------------------------------------------------------
mus=c(0,0,0)
sigmas=c(1,5,10)
x <- seq(0,10,by = 0.5)

list_df <- list()
for (n in seq(1, length(mus), by=1)) {
  list_df[[n]] <- rbind(
    tibble(x = x, n=n, d = distrs[1], 
           y = pfoldnorm(x, mus[n], sigmas[n])),
    tibble(x = x, n=n, d = distrs[2],  
           y = pnorm(x, mus[n], sigmas[n]) - pnorm(-x, mus[n], sigmas[n])),
    tibble(x = x, n=n, d = distrs[3], 
           y = pnorm((x-mus[n])/sigmas[n], 0, 1) - pnorm((-x - mus[n])/sigmas[n], 0, 1)))
}
df <- do.call("rbind", list_df)
df$n <- as.factor(df$n)
df$d <- as.factor(df$d)


p_2d <- ggplot(data=df, aes(x=x,y=y)) +
  geom_line(aes(color= interaction(d,n), linetype=d), lwd=1) +
  ylab("CDF") + xlab("x") +
  theme_classic(base_size=8) + theme(legend.position="none") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_color_manual(values = rep(tri_color,length(mus)))
p_2d

save_plot(paste("figure/F", fig_num, "/F", fig_num, "a integration_cdf_central.tiff", 
                sep = ""), p_2d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# Subfigure B: Compare FN,NN,SN with mu = 1 and several sigmas
#------------------------------------------------------------------------------
mus=c(1,2,4)
sigmas=c(1,1,1)
x <- seq(0,5,by = 0.5)

list_df <- list()
for (n in seq(1, length(mus), by=1)) {
  list_df[[n]] <- rbind(
    tibble(x = x, n=n, d = distrs[1], 
           y = pfoldnorm(x, mus[n], sigmas[n])),
    tibble(x = x, n=n, d = distrs[2],  
           y = pnorm(x, mus[n], sigmas[n]) - pnorm(-x, mus[n], sigmas[n])),
    tibble(x = x, n=n, d = distrs[3], 
           y = pnorm((x-mus[n])/sigmas[n], 0, 1) - pnorm((-x - mus[n])/sigmas[n], 0, 1)))
}
df <- do.call("rbind", list_df)
df$n <- as.factor(df$n)
df$d <- as.factor(df$d)


p_2e <- ggplot(data=df, aes(x=x,y=y)) +
  geom_line(aes(color= interaction(d,n), linetype=d), lwd=1) +
  ylab("CDF") + xlab("x") +
  theme_classic(base_size=8) + theme(legend.position="none") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_color_manual(values = rep(tri_color,length(mus)))
p_2e
save_plot(paste("figure/F", fig_num, "/F", fig_num, "b integration_cdf_central.tiff", 
                sep = ""), p_2e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)





# Generate samples for mu and Sd, calculate all three curves for each, the calculate SSE
n_samples <- 1e4
mus=runif(1e4,-5,5)
sigmas=runif(1e4,.1,5)
x <- seq(0,10,by = 0.5)

df <- tibble( #rmse_FF_vs_FN = rep(0,n_samples),rmse_FF_vs_FS = rep(0,n_samples),
             r_FN_vs_NN= rep(0,n_samples), r_FN_vs_SN= rep(0,n_samples))

for (n in seq(1, n_samples, by=1)) {
  # FN Folded normal
  a = pfoldnorm(x, mus[n], sigmas[n])
  # NN Nonstandard normal
  b = pnorm(x, mus[n], sigmas[n]) - pnorm(-x, mus[n], sigmas[n])
  # SN Standard normal
  c = pnorm( (x-mus[n]) / sigmas[n], 0, 1) - pnorm( (-x-mus[n]) / sigmas[n], 0, 1)
  
  # df$rmse_FF_vs_FN[n] = sqrt(mean((a-b)^2))
  df$r_FN_vs_NN[n] = cor.test(a, b, method = "pearson")$estimate
  
  # df$rmse_FF_vs_FS[n] = sqrt(mean((a-c)^2))
  df$r_FN_vs_SN[n] = cor.test(a, c, method = "pearson")$estimate
}



# QQ plot comparing Folded Normal (FN) with nonstandard normal (NN)
gg <- ggplot(data = tibble(x = a, y = b), aes(x=x, y = y)) +
  geom_point(size=0.5) +
  ylab(expression(CDF~F[NN])) + xlab(expression(CDF~F[FN])) +
  theme_classic(base_size=8) + theme(legend.position="none") +
  geom_abline(slope = 1, intercept = 0, show.legend = NA)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "c Q-Q FN_vs_NN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1.3, dpi = 600)  


# Correlation between FN and NN over large number of QQ plots
gg <- ggplot(data = tibble(y = mean(df$r_FN_vs_NN), sd_y = sd(df$r_FN_vs_NN)),
               aes(x=as.factor(""), y = y)) +
  geom_point(size=0.5) +
  geom_linerange(aes(ymin = y - 1.96*sd_y/sqrt(n_samples), 
                     ymax = y + 1.96*sd_y/sqrt(n_samples))) +
  ylab("Q-Q PCC") + 
  xlab(expression(paste(F[FN], " : ", F[NN], "  ", sep=""))) +
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.text.x=element_blank()) +
  geom_blank(aes(y = 0.985)) +
  geom_blank(aes(y = 1.015))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "d correlation Q-Q FN_vs_NN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 0.75, dpi = 600)  




# QQ plot comparing Folded Normal (FN) with standard normal (SN)
gg <- ggplot(data = tibble(x = a, y = c), aes(x=x, y = y)) +
  geom_point(size=0.5) +
  ylab(expression(CDF~F[SN])) + xlab(expression(CDF~F[FN])) +
  theme_classic(base_size=8) + theme(legend.position="none") +
  geom_abline(slope = 1, intercept =0, show.legend = NA)
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "e Q-Q FN_vs_SN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1.3, dpi = 600)  


# Correlation between FN and NN over large number of QQ plots
gg <- ggplot(data = tibble(y = mean(df$r_FN_vs_SN), sd_y = sd(df$r_FN_vs_SN)),
             aes(x=as.factor(""), y = y)) +
  geom_point(size=0.5) +
  geom_linerange(aes(ymin = y - 1.96*sd_y/sqrt(n_samples), 
                     ymax = y + 1.96*sd_y/sqrt(n_samples))) +
  ylab("Q-Q PCC") + 
  xlab(expression(paste(F[FN], " : ", F[SN], "  ", sep="")))    +
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.text.x=element_blank()) +
  geom_blank(aes(y = 0.985)) +
  geom_blank(aes(y = 1.015))
gg
save_plot(paste("figure/F", fig_num, "/F", fig_num, "f correlation Q-Q FN_vs_SN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = .75, dpi = 600)  