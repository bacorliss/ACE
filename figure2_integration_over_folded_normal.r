

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


fig_basename="f_2"

# dnorm: density function: proability for a specific value of x/zscore of PDF
# pnorm: cum. probability/quantile for a given zscore
# qnorm: return z_score for a given cum. probability/quantile



tri_color <- brewer.pal(n = 3, name = "Set1")


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

save_plot(paste("figure/", fig_basename, "d integration_cdf_central.tiff", 
                sep = ""), p_2d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)


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


save_plot(paste("figure/", fig_basename, "e integration_cdf_central.tiff", 
                sep = ""), p_2e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600)





# Generate samples for mu and Sd, calculate all three curves for each, the calculate SSE
n_samples <- 1e4
mus=runif(1e4,-5,5)
sigmas=runif(1e4,.1,5)
x <- seq(0,10,by = 0.5)

df <- tibble(rmse_FF_FN = rep(0,n_samples),rmse_FF_FS = rep(0,n_samples),
             r_FF_FN= rep(0,n_samples), r_FF_FS= rep(0,n_samples))

for (n in seq(1, n_samples, by=1)) {
  a = pfoldnorm(x, mus[n], sigmas[n])
  b = pnorm(x, mus[n], sigmas[n]) - pnorm(-x, mus[n], sigmas[n])
  c = pnorm( (x-mus[n]) / sigmas[n], 0, 1) - pnorm( (-x-mus[n]) / sigmas[n], 0, 1)
  
  df$rmse_FF_FN[n] = sqrt(mean((a-b)^2))
  df$r_FF_FN[n] = cor.test(a, b, method = "pearson")$estimate
  
  df$rmse_FF_FS[n] = sqrt(mean((a-c)^2))
  df$r_FF_FS[n] = cor.test(a, c, method = "pearson")$estimate
}



p_2e <- ggplot(data=df, aes(y=r_FF_FN)) +
  geom_point() +
  geom_linerange(aes(ymin = cl_fits[1], ymax = cl_fits[2])) +
  ylab("Pearson r") + xlab("") +
  xlab(expression("F[F] : F[N]")) +
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.text.x=element_blank()) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))
p_2e
save_plot(paste("figure/", fig_basename, "r equiv_integration_FF_FN.tiff", 
                sep = ""), p_1f, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1, dpi = 600)  





p_2f <- ggplot(data=df, aes(y=r_FF_FS)) +
  geom_point() +
  geom_linerange(aes(ymin = cl_fits[1], ymax = cl_fits[2])) +
  ylab("Pearson r") + xlab("") +
  xlab(expression("F[F] : F[S]")) +
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.text.x=element_blank()) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))
p_2f
save_plot(paste("figure/", fig_basename, "2f equiv_integration_FF_FS.tiff", 
                sep = ""), p_1f, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1, dpi = 600)  