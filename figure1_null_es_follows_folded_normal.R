

# The magntiude of effect sizefrom normal data follows the folded normal 


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

source("R/stat_helper.r")

# choose colors for plotting
color_pal = brewer.pal(4, "Set1")
fig_num="1"
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), 
           showWarnings = FALSE)


# Parameters
rand_seed <- 0
nsamples = 1e4
equ_font_size=8;
equ_font_size=8
sig_char_size = 3
norm2fnorm_plot_height = 1.75
plot_height = 1.25

# 1A Show examples of histograms of normals transformed into 
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

p_1a <- ggplot(df_1a, aes(x = x)) +
  geom_histogram(aes(fill = group,y=..density..), binwidth = .2, boundary = 0, 
                 position = "identity", alpha = 0.5) + 
  facet_grid(. ~ mu, scales = "free_x", labeller=label_parsed) + 
  ylab("f(x)") +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1a
# Export to TIF
save_plot(paste("figure/F", fig_num, "/F", fig_num, "a example_folded_dists.tiff", sep = ""),
          p_1a, ncol = 1, nrow = 1, base_height = 2,
          base_asp = 3, base_width = 6.5, dpi = 600) 




# B-E) Absolute of normal samples follow folded normal under null hypothesis
#_______________________________________________________________________________
# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data
set.seed(rand_seed)
n_trials <- 1e4
mus <- rep(0, n_trials)#mrunif(1e4, min = -5, max = 5)
n_obs <- 1000

# Generate random samples based on random mu values
# row: ntrials, col: n_obs
x_norm <- t(sapply(mus, function (x) rnorm(n_obs, mean = x, sd = 1),
                   simplify = TRUE))
x_abs_norm <- t(sapply(mus, function (x) abs(rnorm(n_obs, mean = x, sd = 1)),
                       simplify = TRUE))
x_fnorm <- t(sapply(mus, function (x) rfoldnorm(n_obs, mean = x, sd = 1),
                    simplify = TRUE))

# Test if norm is different from abs_norm
p_val_n2abs_n <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_abs_norm[i,]), alternative = "two.sided", 
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)

# Test if norm is different from folded norm
p_val_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                             as.vector(x_fnorm[i,]), alternative = "two.sided", 
                             exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)

# Test if abs_norm is different from folded norm
p_val_abs_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_abs_norm[i,]),
                                 as.vector(x_fnorm[i,]), alternative = "two.sided",
                                 exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
df_1b <- rbind(tibble(d = factor("N : |N|"), x = p_val_n2abs_n),
               tibble(d = factor("N : FN"), x = p_val_n2fn), 
               tibble(d = factor("|N| : FN"), x = p_val_abs_n2fn))

p_1b <- ggplot(df_1b, aes(x=d, y=x)) + 
  geom_violin(scale="width", kernel="gaussian", fill="grey96") +
  xlab("") + ylab("P-Val.") + geom_boxplot(width = 0.1,outlier.size=1) +
  theme_classic(base_size = 8) + theme(legend.position="none", axis.title.x = element_blank())
p_1b
# Export to TIF
save_plot(paste("figure/F", fig_num, "/F", fig_num, "b dist_comparison_violin.tiff", sep = ""),
          p_1b, ncol = 1, nrow = 1, base_height = .8,
          base_asp = 3, base_width = 2.25, dpi = 600)


# Show % of trials that have significant p-values ---------------------------------------------
# Under null distribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_val_n2abs_n<0.05), n_trials, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_val_n2fn<0.05), n_trials, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_val_abs_n2fn<0.05), n_trials, p=0.05, 
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
  geom_blank(aes(y = y_max)) 
p_1c  
save_plot(paste("figure/F", fig_num, "/F", fig_num, "c dist_comparison_facet.tiff",sep = ""), 
          p_1c, ncol = 1, nrow = 1, base_height = .9,
          base_asp = 3, base_width = 1.85, dpi = 600)


# Plot hisogram, of KS values and show that they are derived from a uniform distribution
# Stats
chi_test <- chisq.test(hist(p_val_abs_n2fn, breaks=20)$counts/n_obs, 
                       simulate.p.value = TRUE)
chi_test$p.value
# Plotting
p_1d <- ggplot(tibble(x=p_val_abs_n2fn), aes(x=x)) +
  geom_histogram(aes(y=stat(width*density)), bins=20, fill="grey") +
  theme_classic(base_size=8) + theme(legend.position="none") +
  xlab(expression("P-Values |N| : FN")) + ylab("Probablitly") +
  geom_hline(yintercept=1/20, linetype=1, color="black",alpha=0.5, size=.5)
p_1d             
save_plot(paste("figure/F", fig_num, "/F", fig_num, "d_dist_comparison_hist.tiff", sep = ""), 
          p_1d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)
  


# Plot a single line comparing output between FN and |N|
# Plot difference between P(N) and P(FN) over mu
sort_x_abs_norm = sort(x_abs_norm[1,])
sort_x_fnorm = sort(x_fnorm[1,])

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
fit <- unname(lm(formula = sort_x_abs_norm ~ 0 +sort_x_fnorm)$coefficients[1])

# Define data for plotting
df_1e = tibble(x = sort_x_fnorm, y = sort_x_abs_norm)
p_1e <- ggplot(df_1e, aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab(expression("Sorted Sample FN")) + ylab("Sorted Sample |N|   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1e
save_plot(paste("figure/F", fig_num, "/F", fig_num, "e dist_comp_qq_plot.tiff",sep = ""), 
          p_1e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)

# Fit a line for each of the simulations
lin_model <- function(x,y) lm(formula = y ~ 0 + x)$coefficients[1];
fits <- unname(sapply(1:dim(x_fnorm)[1], 
                      function(i) lin_model(sort(x_fnorm[i,]), sort(x_abs_norm[i,])),
                      simplify = TRUE ))
cl_fits <- quantile(fits, c(.025, .975)) 
df_fits <- tibble(x=as.factor("|N| : FN"),y=mean(fits),y_min = 0.99,y_max = 1.1)

p_1f <- ggplot(data=df_fits, aes(x=x, y=y)) +
  geom_point() +
  #geom_hline(yintercept = 1, linetype = 2, color = "red", alpha = .2) +
  geom_linerange(aes(ymin = cl_fits[1], ymax = cl_fits[2])) +
  ylab("Q-Q Slopes") + xlab("") +
  #xlab(expression("|N| : FN")) +
  theme_classic(base_size=8) + theme(legend.position="none") +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))
p_1f
save_plot(paste("figure/F", fig_num, "/F", fig_num, "f dist_comparison_qq_slope.tiff", 
                sep = ""), p_1f, ncol = 1, nrow = 1, base_height = 1.45,
          base_asp = 3, base_width = 1, dpi = 600)





# G-F) Show that absolute of normal samples follow folded normal under alternate hypothesis
#______________________________________________________________________________________
# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data
set.seed(rand_seed)
n_trials <- 1e4
mus <- runif(1e4, min = -5, max = 5)
sds <-  runif(1e4, min = .1, max = 5)
n_obs <- 1000

# Generate random samples based on random mu values
x_norm <- t(sapply(mus, function (x) rnorm(n_obs, mean = x, sd = 1), 
                   simplify = TRUE))
x_abs_norm <- t(sapply(mus, function (x) abs(rnorm(n_obs, mean = x, sd = 1)), 
                       simplify = TRUE))
x_fnorm <- t(sapply(mus, function (x) rfoldnorm(n_obs, mean = x, sd = 1), 
                    simplify = TRUE))

# Test if norm is different from abs_norm
p_val_n2abs_n <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_abs_norm[i,]), alternative = "two.sided",
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
# Test if norm is different from folded norm
p_val_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_fnorm[i,]), alternative = "two.sided", 
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
# Test if abs_norm is different from folded norm
p_val_abs_n2fn <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_abs_norm[i,]),
                                as.vector(x_fnorm[i,]), alternative = "two.sided", 
                                exact = TRUE, tol=1e-8, simulate.p.value=FALSE)$p)
df_1g <- rbind(tibble(d = factor("N : |N|"), x = p_val_n2abs_n),
               tibble(d = factor("N : FN"), x = p_val_n2fn), 
               tibble(d = factor("|N| : FN"), x = p_val_abs_n2fn))

p_1g <- ggplot(df_1g, aes(x=d, y=x)) + 
  geom_violin(scale="width", kernel="gaussian", fill="grey96") +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  xlab("") + ylab("P-Val.") + 
  theme_classic(base_size=8) + theme(legend.position="none",
                                     axis.title.x=element_blank())
p_1g
# Export to TIF
save_plot(paste("figure/F", fig_num, "/F", fig_num, "g dist_comparison_violin.tiff", sep = ""), 
          p_1g, ncol = 1, nrow = 1, base_height = .75,
          base_asp = 3, base_width = 2.25, dpi = 600)


# Show % of trials that have significant p-values ---------------------------------------------
# Under null disctribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_val_n2abs_n<0.05), n_trials, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_val_n2fn<0.05), n_trials, p=0.05, 
                      conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_val_abs_n2fn<0.05), n_trials, p=0.05, 
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
  geom_blank(aes(y = y_max)) 
p_1h  
save_plot(paste("figure/F", fig_num, "/F", fig_num, "h dist_comparison_facet.tiff",sep = ""), 
          p_1h, ncol = 1, nrow = 1, base_height = 1,
          base_asp = 3, base_width = 1.85, dpi = 600)


# Plot hisogram, of KS values and show that they are derived from a uniform distribution
# Stats
chi_test <- chisq.test(hist(p_val_abs_n2fn, breaks=20)$counts/n_obs, 
                       simulate.p.value = TRUE)
chi_test$p.value
# Plotting
p_1i <- ggplot(tibble(x=p_val_abs_n2fn), aes(x=x)) +
  geom_histogram(aes(y=stat(width*density)), bins=20, fill="grey") +
  theme_classic(base_size=8) + theme(legend.position="none") +
  xlab("P-Values |N| : FN") + ylab("Probablitly") +
  geom_hline(yintercept=1/20, linetype=1, color="black",alpha=0.5, size=.5)
p_1i            
save_plot(paste("figure/F", fig_num, "/F", fig_num, "i_dist_comparison_hist.tiff", sep = ""), 
          p_1i, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)



# Plot a single line comparing output between FN and |N|
# Plot difference between P(N) and P(F[N]) over mu
sort_x_abs_norm = sort(x_abs_norm[1,])
sort_x_fnorm = sort(x_fnorm[1,])

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
fit <- unname(lm(formula = sort_x_abs_norm ~ 0 +sort_x_fnorm)$coefficients[1])

# Define data for plotting
df_1j = tibble(x=sort_x_fnorm, y=sort_x_abs_norm)
p_1j <- ggplot(df_1j, aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab("Sorted Sample FN") + ylab("Sorted Sample |N|   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_classic(base_size=8) + theme(legend.position="none")
p_1j
save_plot(paste("figure/F", fig_num, "/F", fig_num, "j dist_comp_qq_plot.tiff",sep = ""), 
          p_1j, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600)

# Fit a line for each of the simulations
lin_model <- function(x,y) lm(formula = y ~ 0 + x)$coefficients[1];
fits <- unname(sapply(1:dim(x_fnorm)[1], 
                      function(i) lin_model(sort(x_fnorm[i,]), sort(x_abs_norm[i,])),
                      simplify = TRUE ))
cl_fits <- quantile(fits, c(.025, .975)) 
df_fits <- tibble(x=as.factor("|N|:FN"),y=mean(fits),y_min = 0.99,y_max = 1.1)

p_1k <- ggplot(data=df_fits, aes(x=x,y=y)) +
  geom_point() +
  #geom_hline(yintercept = 1, linetype = 2, color = "red", alpha = .2) +
  geom_linerange(aes(ymin = cl_fits[1], ymax = cl_fits[2])) +
  ylab("Q-Q Slopes") + xlab("") +
  theme_classic(base_size=8) + theme(legend.position="none") +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) 
p_1k
save_plot(paste("figure/F", fig_num, "/F", fig_num, "k dist_comparison_qq_slope.tiff", 
                sep = ""), p_1k, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1, dpi = 600)

