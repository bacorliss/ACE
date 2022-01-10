
# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)


# Load packages
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(broom)
p_load(gridExtra)
p_load(grid)
p_load(rlang)
p_load(colorspace)
p_load(VGAM)
p_load(boot)
p_load(dplyr)
p_load(cowplot)
base_dir = "mdm_t"

fig_num = "2"
dir.create(file.path(getwd(), paste(base_dir,"/figure/F",fig_num,sep="")), 
           showWarnings = FALSE, recursive = TRUE)


# dnorm: density function: proability for a specific value of x/zscore of PDF
# pnorm: cum. probability/quantile for a given zscore
# qnorm: return z_score for a given cum. probability/quantile


# User functions
#------------------------------------------------------------------------------
# CDF for folded normal, three different functions
# Folded normal
fun_pfoldnorm <- function (x, mu, sigma) {
  p = pnorm(x,mean=mu, sd=sigma) - pnorm(0,mean=mu, sd=sigma) +
    pnorm(0,mean=mu, sd=sigma) - pnorm(-x,mean=mu, sd=sigma)
}
# Nonstandard Normal
fun_pnstandnorm <- function (x, mu, sigma) {pnorm(x, mean = mu, sd = sigma) -
    pnorm(-x, mu, sigma)}
# Standard Normal
fun_pstandnorm <- function (x, mu, sigma) {pnorm((x-mu)/sigma, mean = 0, sd = 1) -
    pnorm((-x - mu)/sigma, mean = 0, sd = 1)}


# Script parameters
tri_color <- brewer.pal(n = 3, name = "Set1")
distrs=c("FN", "NN","SN")
rand.seed <- 0



# Subfigure A: Compare FN,NN,SN with mu = 1 and swept sigma
#------------------------------------------------------------------------------
mus = c(0,0,0)
sigmas = c(1,5,10)
x <- seq(0,10,by = 0.5)

list_df <- list()
for (n in seq(1, length(mus), by=1)) {
  list_df[[n]] <- rbind(
    tibble(x = x, n=n, d = distrs[1], 
           y = fun_pfoldnorm(x, mu = mus[n], sigma = sigmas[n])),
    tibble(x = x, n = n, d = distrs[2],  
           y = fun_pnstandnorm(x, mu = mus[n], sigma = sigmas[n])),
    tibble(x = x, n=n, d = distrs[3], 
           y = fun_pstandnorm(x, mu = mus[n], sigma = sigmas[n])))
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

save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "a integration_cdf_central.tiff", 
                sep = ""), p_2d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


# Subfigure B: Compare FN,NN,SN with mu = 1 and several sigmas
#------------------------------------------------------------------------------
mus = c(1,2,4)
sigmas = c(1,1,1)
x <- seq(0,5,by = 0.5)

list_df <- list()
for (n in seq(1, length(mus), by=1)) {
  list_df[[n]] <- rbind(
    tibble(x = x, n=n, d = distrs[1], 
           y = fun_pfoldnorm(x, mu = mus[n], sigma = sigmas[n])),
    tibble(x = x, n = n, d = distrs[2],  
           y = fun_pnstandnorm(x, mu = mus[n], sigma = sigmas[n])),
    tibble(x = x, n=n, d = distrs[3], 
           y = fun_pstandnorm(x, mu = mus[n], sigma = sigmas[n])))
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
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "b integration_cdf_central.tiff", 
                sep = ""), p_2e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)





# Generate samples for mu and Sd, calculate all three curves for each, the calculate SSE
set.seed(rand.seed)
n_samples <- 1e4
mus = runif(n_samples,-5,5)
sigmas = runif(n_samples,.1,5)
x <- seq(0,10,by = 0.5)

df <- tibble( #rmse_FF_vs_FN = rep(0,n_samples),rmse_FF_vs_FS = rep(0,n_samples),
             r_FN_vs_NN= rep(0,n_samples), r_FN_vs_SN= rep(0,n_samples),
             a_FN_vs_NN= rep(0,n_samples), a_FN_vs_SN= rep(0,n_samples))

for (n in seq(1, n_samples, by=1)) {
  # FN Folded normal
  F_FN = pfoldnorm(x, mus[n], sigmas[n])
  # NN Nonstandard normal
  F_NN = pnorm(x, mus[n], sigmas[n]) - pnorm(-x, mus[n], sigmas[n])
  # SN Standard normal
  F_SN = pnorm( (x-mus[n]) / sigmas[n], 0, 1) - pnorm( (-x-mus[n]) / sigmas[n], 0, 1)
  
  # df$rmse_FF_vs_FN[n] = sqrt(mean((a-b)^2))
  # df$rmse_FF_vs_FS[n] = sqrt(mean((a-c)^2))
  df$r_FN_vs_NN[n] = cor.test(F_FN, F_NN, method = "pearson")$estimate
  df$r_FN_vs_SN[n] = cor.test(F_FN, F_SN, method = "pearson")$estimate
  
  lm_FN_vs_NN <- lm(F_FN~0+F_NN)
  df$a_FN_vs_NN[n] = lm_FN_vs_NN$coefficients[1]
  
  lm_FN_vs_SN <- lm(F_FN~0+F_SN)
  df$a_FN_vs_SN[n] = lm_FN_vs_SN$coefficients[1]

}



# QQ plot comparing Folded Normal (FN) with nonstandard normal (NN)
gg <- ggplot(data = tibble(x = F_FN, y = F_NN), aes(x=x, y = y)) +
  geom_point(size=0.5) +
  ylab(expression(Prob.~F[NN])) + xlab(expression(Prob.~F[FN])) +
  theme_classic(base_size=8) + theme(legend.position="none") +
  geom_abline(slope = 1, intercept = 0, show.legend = NA)
gg
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "c Q-Q FN_vs_NN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1.3, dpi = 600)  


# Correlation between FN and NN over large number of QQ plots
gg <- ggplot(data = tibble(y = mean(df$r_FN_vs_NN), sd_y = sd(df$r_FN_vs_NN)),
               aes(x=as.factor(""), y = y)) +
  geom_point(size=0.5) +
  geom_linerange(aes(ymin = y - 1.96*sd_y/sqrt(n_samples), 
                     ymax = y + 1.96*sd_y/sqrt(n_samples))) +
  ylab("Slope") + 
  xlab(expression(paste(F[FN], " : ", F[NN], "    ", sep=""))) +
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.text.x=element_blank())
gg
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "d correlation Q-Q FN_vs_NN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 0.75, dpi = 600)  




# QQ plot comparing Folded Normal (FN) with standard normal (SN)
gg <- ggplot(data = tibble(x = F_FN, y = F_SN), aes(x=x, y = y)) +
  geom_point(size=0.5) +
  ylab(expression(Prob.~F[SN])) + xlab(expression(Prob.~F[FN])) +
  theme_classic(base_size=8) + theme(legend.position="none") +
  geom_abline(slope = 1, intercept =0, show.legend = NA)
gg
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "e Q-Q FN_vs_SN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1.3, dpi = 600)  


# Correlation between FN and NN over large number of QQ plots
gg <- ggplot(data = tibble(y = mean(df$a_FN_vs_SN), sd_y = sd(df$a_FN_vs_SN)),
             aes(x=as.factor(""), y = y)) +
  geom_point(size=0.5) +
  geom_linerange(aes(ymin = y - 1.96*sd_y/sqrt(n_samples), 
                     ymax = y + 1.96*sd_y/sqrt(n_samples))) +
  ylab("Slope") + 
  xlab(expression(paste(F[FN], " : ", F[SN], "    ", sep="")))    +
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.text.x=element_blank())
  # geom_blank(aes(y = 1+1e-17)) +
  # geom_blank(aes(y = 1-1e-17))
gg
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "f correlation Q-Q FN_vs_SN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = .75, dpi = 600)  










# Calculate MDM with three different integration methods
#-------------------------------------------------------------------------------
n_samples <- 1e3
mus=runif(n_samples,-5,5)
sigmas=runif(n_samples,.1,5)
conf.level = 0.95


mdm_fn   <- rep(0,length(mus))
mdm_nn <- rep(0,length(mus))
mdm_sn <- rep(0,length(mus))
for (n in seq_along(mus)) {
  mdm_fn [n] <- uniroot(function(x) fun_pfoldnorm(x, mus[n],sigmas[n])-0.95, 
                              c(abs(mus[n]), abs(mus[n]) + 6*sigmas[n]),  
                              tol = .Machine$double.eps)$root
  mdm_nn [n] <- uniroot(function(x) fun_pnstandnorm(x, mus[n],sigmas[n])-0.95, 
                                c(abs(mus[n]), abs(mus[n]) + 6*sigmas[n]),  
                                tol = .Machine$double.eps)$root
  mdm_sn [n] <- uniroot(function(x) fun_pstandnorm(x, mus[n],sigmas[n])-0.95, 
                               c(abs(mus[n]), abs(mus[n]) + 6*sigmas[n]),  
                               tol = .Machine$double.eps)$root
}

# Compare FN vs NN
df_fn_vs_nn <- tibble(x = (mdm_fn + mdm_nn)/2, y = (mdm_fn - mdm_nn) / 
                        (mdm_fn + mdm_nn)/2)
gg <- ggplot(df_fn_vs_nn, aes(x=x,y=y)) +
  geom_hline(yintercept = 1.96*sd(df_fn_vs_nn$y), color = "red", linetype="dashed", size=0.25) +
  geom_hline(yintercept = -1.96*sd(df_fn_vs_nn$y), color = "red", linetype="dashed", size=0.25) +
  geom_hline(yintercept = 0, color="blue", size=0.25)+
  geom_point(size=0.1) +
  xlab('Mean') + 
  ylab('Rel. Diff') +
  theme_classic(base_size=8) +
  geom_blank(aes(y = -0.6E-15)) +
  geom_blank(aes(y = .6E-15))
gg
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "g_bland_altman MDM_FN_vs_NN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 2, dpi = 600)  
# Compare FN vs SN
df_fn_vs_sn <- tibble(x = (mdm_fn + mdm_sn)/2, y = (mdm_fn - mdm_sn) / 
                        (mdm_fn + mdm_sn)/2)
gg <- ggplot(df_fn_vs_nn, aes(x=x,y=y)) +
  geom_hline(yintercept = 1.96*sd(df_fn_vs_sn$y), color = "red", linetype="dashed", size=0.25) +
  geom_hline(yintercept = -1.96*sd(df_fn_vs_sn$y), color = "red", linetype="dashed", size=0.25) +
  geom_hline(yintercept = 0, color="blue", size=0.25)+
  geom_point(size=0.1) +
  xlab('Mean') + 
  ylab('Rel. Diff') +
  theme_classic(base_size=8) +
  geom_blank(aes(y = -0.6E-15)) +
  geom_blank(aes(y = 0.6E-15))
gg
save_plot(paste(base_dir,"/figure/F", fig_num, "/F", fig_num, "g_bland_altman MDM_FN_vs_SN.tiff", 
                sep = ""), gg, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 2, dpi = 600) 







