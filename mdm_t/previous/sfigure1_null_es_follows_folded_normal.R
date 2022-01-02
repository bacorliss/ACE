

#' Empirical evidence that that agreement of means follows a folded normal


# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}
library(pacman)
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
p_load(VGAM)
# source("R/stat_helper.r")
base_dir = "mdm_t"


# Figure parameters
#-------------------------------------------------------------------------------
color_pal = brewer.pal(4, "Set1")
fig_num="1"
dir.create(file.path(getwd(), paste(base_dir, "/figure/SF",fig_num,sep="")), 
           recursive = TRUE, showWarnings = FALSE)


# Parameters
rand_seed <- 0
nsamples = 1e4
equ_font_size=8;
equ_font_size=8
sig_char_size = 3
norm2fnorm_plot_height = 1.75
plot_height = 1.25

# Figure 1A: Example of normal distribution before and after absolute folding
# at various mu values
set.seed(rand_seed)
sample_1a = rnorm(nsamples, mean = 0, sd = 1)
df_1a = rbind(tibble(group = 2, mu = "mu==0", x = sample_1a), 
                tibble(group = 1, mu = "mu==0", x = abs(sample_1a)),
              tibble(group = 2, mu = "mu==0.5", x = sample_1a + 0.5), 
                tibble(group = 1, mu = "mu==0.5", x = abs(sample_1a + 0.5)),
              tibble(group = 2, mu = "mu==1", x = sample_1a + 1), 
                tibble(group = 1, mu = "mu==1", x = abs(sample_1a + 1)),
              tibble(group = 2, mu = "mu==1.5", x = sample_1a + 1.5), 
                tibble(group = 1, mu = "mu==1.5", x = abs(sample_1a + 1.5)),
              tibble(group = 2, mu = "mu==2", x = sample_1a + 2), 
                tibble(group = 1, mu = "mu==2", x = abs(sample_1a + 2)),
              #tibble(group = 2, mu = "mu==2.5", x = sample_1a + 2.5), 
              #  tibble(group = 1, mu = "mu==2.5", x = abs(sample_1a + 2.5)),
              tibble(group = 2, mu = "mu==3", x = sample_1a + 3), 
                tibble(group = 1, mu = "mu==3", x = abs(sample_1a + 3))
              )
df_1a$group <- as.factor(df_1a$group)
df_1a$mu    <- as.factor(df_1a$mu)

df_1a %>%
  group_by(group,mu) %>%
  summarise(mean = mean(x), n = n())

df_summary <- df_1a %>%
  group_by(mu) %>%
  summarize(ks_p_value = min(ks.test(x[group == 1], x[group == 2])$p.value*6,1))




p_1a <- ggplot(df_1a, aes(x = x)) +
  geom_histogram(aes(fill = group,y=..density..), binwidth = .2, boundary = 0, 
                 position = "identity", alpha = 0.5) + 
  facet_grid(. ~ mu, scales = "free_x", labeller=label_parsed) + 
  # geom_text(label = df_summary$ks_p_value, aes(y=0.07)) + 
  ylab("f(x)") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1a
# Export to TIF 
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "a example_folded_dists.tiff", sep = ""),
          p_1a, ncol = 1, nrow = 1, base_height = 1.2,
          base_asp = 3, base_width = 6.5, dpi = 600) 




# Figure B:  Absolute of normal samples follow folded normal under null hypothesis
#_______________________________________________________________________________
# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data
set.seed(rand_seed)
n_sims <- 1e4
mus <- rep(0, n_sims)
n_obs <- 1000

# Generate random samples based on random mu values
# row: ntrials, col: n_obs
x_norm <- t(sapply(mus, function (x) rnorm(n_obs, mean = x, sd = 1),
                   simplify = TRUE))
x_absnorm <- t(sapply(mus, function (x) abs(rnorm(n_obs, mean = x, sd = 1)),
                       simplify = TRUE))
x_fnorm <- t(sapply(mus, function (x) rfoldnorm(n_obs, mean = x, sd = 1),
                    simplify = TRUE))

# Test if norm is different from absnorm
p_val_n2abs_n <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_absnorm[i,]), alternative = "two.sided", 
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)

# Test if norm is different from folded norm
p_val_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                             as.vector(x_fnorm[i,]), alternative = "two.sided", 
                             exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)

# Test if absnorm is different from folded norm
p_val_abs_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_absnorm[i,]),
                                 as.vector(x_fnorm[i,]), alternative = "two.sided",
                                 exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
df_1b <- rbind(tibble(d = factor("N : |N|"), x = p_val_n2abs_n),
               tibble(d = factor("N : FN"), x = p_val_n2fn), 
               tibble(d = factor("|N| : FN"), x = p_val_abs_n2fn))

# Test if each p-value is from normal distribution
unif_p_vals = c(ifelse(length(unique(p_val_n2abs_n))==1, "p < 1.2e-6",
                       chisq.test(hist(p_val_n2abs_n, breaks=20)$counts/n_obs, 
                                  simulate.p.value = FALSE)$p.value),
                ifelse(length(unique(p_val_n2fn))==1, "p < 1.2e-6",
                       chisq.test(hist(p_val_n2fn, breaks=20)$counts/n_obs, 
                           simulate.p.value = FALSE)$p.value),
                ifelse(length(unique(p_val_abs_n2fn))==1, "p < 1.e-6",
                       sprintf("p = %.3f", chisq.test(hist(p_val_abs_n2fn, breaks=20)$counts/n_obs, 
                           simulate.p.value = FALSE)$p.value)))



p_1b <- ggplot(df_1b, aes(x=d, y=x)) + 
  geom_violin(scale="width", kernel="gaussian", fill="grey96") +
  geom_text(data = tibble(d = levels(df_1b$d),y = rep(1.15,3)), aes(x=d,y=y, label = unif_p_vals), vjust = 0.5, size=1.8) +
  geom_hline(yintercept=1, size=0.2,alpha = 0.2)+
  xlab("") + ylab("KS P-Val.") + geom_boxplot(width = 0.1,outlier.size=1) +
  theme_classic(base_size = 8) + theme(legend.position="none", axis.title.x = element_blank())
p_1b
# Export to TIF
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "b dist_comparison_violin.tiff", sep = ""),
          p_1b, ncol = 1, nrow = 1, base_height = .8,
          base_asp = 2, base_width = 3, dpi = 600)


# Show % of trials that have significant p-values ---------------------------------------------
# Under null distribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_val_n2abs_n<0.05), n_sims, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_val_n2fn<0.05), n_sims, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_val_abs_n2fn<0.05), n_sims, p=0.05, 
                       conf.level=1-0.05/3, correct = FALSE)

comp_levels <- c("N : |N|","N : FN","|N| : FN")
df_1c <- tibble(d = factor(comp_levels,levels = comp_levels, ordered = TRUE), 
       estimate = c(bin_n_an$estimate,bin_n_fn$estimate,bin_an_fn$estimate),
       lcl = c(bin_n_an$conf.int[1],bin_n_fn$conf.int[1],bin_an_fn$conf.int[1]),
       ucl = c(bin_n_an$conf.int[2],bin_n_fn$conf.int[2],bin_an_fn$conf.int[2]),
       p_value = c(bin_n_an$p.value,bin_n_fn$p.value,bin_an_fn$p.value),
       y_min = c(.997,.997,0.04),y_max = c(1, 1, 0.06),
       y_horz_line=c(1,1,0.05))


p_1c <- ggplot(df_1c, aes(x=d, y=estimate)) + 
  geom_point(size=1) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1) +
  facet_wrap(~d,scales="free", strip.position = 'bottom') +
  xlab("") +  ylab("Trials p<\u03B1 ") +
  theme_classic(base_size=8) + 
  theme(legend.position="none",strip.text.x = element_blank(), 
        axis.ticks.x=element_blank()) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) +
scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
p_1c  
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "c dist_comparison_facet.tiff",sep = ""), 
          p_1c, ncol = 1, nrow = 1, base_height = .9,
          base_asp = 3, base_width = 2.6, dpi = 600)


# Plot a single line comparing output between FN and |N|
# Plot difference between P(N) and P(FN) over mu
sort_x_norm <- sort(x_norm[1,])
sort_x_absnorm = sort(x_absnorm[1,])
sort_x_fnorm = sort(x_fnorm[1,])


q_sort_x_absnorm <- cumsum(sort_x_absnorm)/sum(sort_x_absnorm)
q_sort_x_fnorm    <- cumsum(sort_x_fnorm)/sum(sort_x_fnorm)


## TODO Change QQ plot to compare sample to theoretical quantile, not working yet
# abs_sort_x_norm <- abs(sort_x_norm)
# z_sort_x_norm <- scale(sort_x_norm, center = TRUE, scale = TRUE)
# q_abs_sort_x_norm <- cumsum(abs_sort_x_norm)/sum(abs_sort_x_norm)
# qt_fnorm <-qfoldnorm(abs(z_sort_x_norm))
# plot(qt_fnorm,q_abs_sort_x_norm)

# Calculate theoretical quantiles of x_absnorm  
# qfoldnorm(scale(sort_x_absnorm, center = TRUE, scale = TRUE))

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
# fit <- unname(lm(formula = sort_x_absnorm ~ 0 +sort_x_fnorm)$coefficients[1])




# Define data for plotting
p_1e <- ggplot(tibble(x = q_sort_x_fnorm, y = q_sort_x_absnorm), aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab(expression("Sample Quantiles FN")) + ylab("Sample Quantile |N|   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1e
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "e dist_comp_qq_plot.tiff",sep = ""), 
          p_1e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)

# Fit a line for each of the simulations
lin_model <- function(x,y) lm(formula = y ~ 0 + x)$coefficients[1];
fits <- unname(sapply(1:dim(x_fnorm)[1], 
                      function(i) lin_model(
                        cumsum(sort(x_fnorm[i,]))   /sum(x_fnorm[i,]), 
                        cumsum(sort(x_absnorm[i,]))/sum(x_fnorm[i,]) ),
                      simplify = TRUE ))
ttest_fits <- t.test(fits, mu = 1)


p_1f <- ggplot(data=tibble(x=as.factor("|N| : FN"),y = mean(fits)), aes(x=x, y=y)) +
  geom_point() +
  geom_linerange(aes(ymin =  ttest_fits$conf.int[1], ymax =  ttest_fits$conf.int[2])) +
  ylab("Q-Q Slopes") + xlab("") +
  #xlab(expression("|N| : FN")) +
  theme_classic(base_size=8) + theme(legend.position="none")
  # geom_blank(aes(y = 0.999)) +
  # geom_blank(aes(y = 1.001))
p_1f
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "f dist_comparison_qq_slope.tiff", 
                sep = ""), p_1f, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1, dpi = 600)





# G-F) Show that absolute of normal samples follow folded normal under alternate hypothesis
#______________________________________________________________________________________
# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data
set.seed(rand_seed)
mus <- runif(1e4, min = -5, max = 5)
sds <-  runif(1e4, min = .1, max = 5)
n_obs <- 1000

# Generate random samples based on random mu values
x_norm <- t(sapply(mus, function (x) rnorm(n_obs, mean = x, sd = 1), 
                   simplify = TRUE))
x_absnorm <- t(sapply(mus, function (x) abs(rnorm(n_obs, mean = x, sd = 1)), 
                       simplify = TRUE))
x_fnorm <- t(sapply(mus, function (x) rfoldnorm(n_obs, mean = x, sd = 1), 
                    simplify = TRUE))

# Test if norm is different from absnorm
p_val_n2abs_n <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_absnorm[i,]), alternative = "two.sided",
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
# Test if norm is different from folded norm
p_val_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_fnorm[i,]), alternative = "two.sided", 
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
# Test if absnorm is different from folded norm
p_val_abs_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_absnorm[i,]),
                                as.vector(x_fnorm[i,]), alternative = "two.sided", 
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
df_1g <- rbind(tibble(d = factor("N : |N|"), x = p_val_n2abs_n),
               tibble(d = factor("N : FN"), x = p_val_n2fn), 
               tibble(d = factor("|N| : FN"), x = p_val_abs_n2fn))

# Test if each p-value is from normal distribution
unif_p_vals = c(ifelse(length(unique(p_val_n2abs_n))==1, "p < 1.2e-6",
                       sprintf("p = %.2e", chisq.test(hist(p_val_n2abs_n, breaks=20)$counts/n_obs, 
                                  simulate.p.value = FALSE)$p.value)),
                ifelse(length(unique(p_val_n2fn))==1, "p < 1.2e-6",
                       sprintf("p = %.2e", chisq.test(hist(p_val_n2fn, breaks=20)$counts/n_obs, 
                                  simulate.p.value = FALSE)$p.value)),
                ifelse(length(unique(p_val_abs_n2fn))==1, "p < 1.e-6",
                       sprintf("p = %.2e", chisq.test(hist(p_val_abs_n2fn, breaks=20)$counts/n_obs, 
                                                      simulate.p.value = FALSE)$p.value)))


p_1g <- ggplot(df_1g, aes(x=d, y=x)) + 
  geom_violin(scale="width", kernel="gaussian", fill="grey96") +
  geom_text(data = tibble(d = levels(df_1b$d),y = rep(1.15,3)), aes(x=d,y=y, label = unif_p_vals), vjust = 0.5, size=1.8) +
  geom_hline(yintercept=1, size=0.2,alpha = 0.2)+
  geom_boxplot(width=0.1, outlier.shape=NA) +
  xlab("") + ylab("P-Val.") + 
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.title.x=element_blank())
p_1g
# Export to TIF
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "g dist_comparison_violin.tiff", sep = ""), 
          p_1g, ncol = 1, nrow = 1, base_height = .75,
          base_asp = 3, base_width = 3, dpi = 600)


# Show % of trials that have significant p-values ---------------------------------------------
# Under null disctribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_val_n2abs_n<0.05), n_sims, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_val_n2fn<0.05), n_sims, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_val_abs_n2fn<0.05), n_sims, p=0.05, 
                       conf.level=1-0.05/3, correct = FALSE)

comp_levels <- c("N : |N|","N : FN","|N| : FN")
df_1h <- tibble(d = factor(comp_levels,levels = comp_levels, ordered = TRUE), 
                estimate = c(bin_n_an$estimate,bin_n_fn$estimate,bin_an_fn$estimate),
                lcl = c(bin_n_an$conf.int[1],bin_n_fn$conf.int[1],bin_an_fn$conf.int[1]),
                ucl = c(bin_n_an$conf.int[2],bin_n_fn$conf.int[2],bin_an_fn$conf.int[2]),
                p_value = c(bin_n_an$p.value,bin_n_fn$p.value,bin_an_fn$p.value),
                y_min = c(.997,.997,0.04),y_max = c(1, 1, 0.06),
                y_horz_line=c(1,1,0.05))


p_1h <- ggplot(df_1h, aes(x=d, y=estimate)) + 
  geom_point(size=1) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1) +
  facet_wrap(~d,scales="free", strip.position = 'bottom') +
  xlab("") +  ylab("Trials p<\u03B1 ") +
  theme_classic(base_size=8) + theme(legend.position="none",strip.text.x = element_blank(),
                                     axis.ticks.x=element_blank()) + #axis.text.x=element_blank()
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
p_1h  
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "h dist_comparison_facet.tiff",sep = ""), 
          p_1h, ncol = 1, nrow = 1, base_height = 1,
          base_asp = 3, base_width = 2.5, dpi = 600)





# Plot a single line comparing output between FN and |N|
# Plot difference between P(N) and P(F[N]) over mu
sort_x_absnorm = sort(x_absnorm[1,])
sort_x_fnorm = sort(x_fnorm[1,])

q_sort_x_absnorm <- cumsum(sort_x_absnorm)/sum(sort_x_absnorm)
q_sort_x_fnorm    <- cumsum(sort_x_fnorm)/sum(sort_x_fnorm)

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
fit <- unname(lm(formula = q_sort_x_absnorm ~ 0 +q_sort_x_fnorm)$coefficients[1])

# Define data for plotting
df_1j = tibble(x=q_sort_x_fnorm, y=q_sort_x_absnorm)
p_1j <- ggplot(df_1j, aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab("Sample Quantile FN") + ylab("Sorted Sample |N|   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1j
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "j dist_comp_qq_plot.tiff",sep = ""), 
          p_1j, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)

# Fit a line for each of the simulations
lin_model <- function(x,y) lm(formula = y ~ 0 + x)$coefficients[1];
fits <- unname(sapply(1:dim(x_fnorm)[1], 
                      function(i) lin_model(
                        cumsum(sort(x_fnorm[i,]))   /sum(x_fnorm[i,]), 
                        cumsum(sort(x_absnorm[i,]))/sum(x_fnorm[i,])),
                      simplify = TRUE ))
ttest_fits <- t.test(fits, mu = 1)

# cl_fits <- quantile(fits, c(.025, .975)) 
p_1k <- ggplot(data=tibble(x=as.factor("|N| : FN"),y = mean(fits)), aes(x=x, y=y)) +
  geom_point() +
  geom_linerange(aes(ymin = ttest_fits$conf.int[1], ymax = ttest_fits$conf.int[2])) +
  ylab("Q-Q Slopes") + xlab("") +
  theme_classic(base_size=8) + theme(legend.position="none")
  # geom_blank(aes(y = 0.999)) +
  # geom_blank(aes(y = 1.001)) 
p_1k
save_plot(paste(base_dir, "/figure/SF", fig_num, "/F", fig_num, "k dist_comparison_qq_slope.tiff", 
                sep = ""), p_1k, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1, dpi = 600)

