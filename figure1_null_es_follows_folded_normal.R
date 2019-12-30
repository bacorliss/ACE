

# The magntiude of effect sizefrom normal data follows the folded normal 
# distribution and converges to the normal distribution as mu>>sigma




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

# choose colors for plotting
color_pal = brewer.pal(4, "Set1")


# Parameters
set.seed(0)
dist_text_size=3;
equ_font_size=8;

# Initial parameters
nsamples = 1e6
ngroups = 2
set.seed(0)
dist_text_size=3;
equ_font_size=8;


theme_small_classic <- function () 
{ 
  theme_classic(base_size=8)
}




# 1A Show examples of histograms of normals transformed into 
sample_1a = rnorm(nsamples, mean = 0, sd = 1)
df_1a = rbind(tibble(group = 2, mu = "mu==0", x = sample_1a    ), tibble(group = 1, mu = "mu==0", x = abs(sample_1a)    ),
              tibble(group = 2, mu = "mu==1", x = sample_1a + 1), tibble(group = 1, mu = "mu==1", x = abs(sample_1a + 1)),
              tibble(group = 2, mu = "mu==1.5", x = sample_1a + 1.5), tibble(group = 1, mu = "mu==1.5", x = abs(sample_1a + 1.5)),
              tibble(group = 2, mu = "mu==3", x = sample_1a + 3), tibble(group = 1, mu = "mu==3", x = abs(sample_1a + 3))
              )
df_1a$group <- as.factor(df_1a$group)
df_1a$mu    <- as.factor(df_1a$mu)

df_1a %>%
  group_by(group,mu) %>%
  summarise(mean = mean(x), n = n())

p_1a <- ggplot(df_1a, aes(x = x)) +
  geom_histogram(aes(fill = group,y=..density..), binwidth = .2, boundary = 0, position = "identity", alpha = 0.5) + 
  facet_grid(. ~ mu, scales = "free_x", labeller=label_parsed) + ylab("Probability") +
  theme_small_classic()
p_1a
# Export to TIF
save_plot("figure/figure_1a_example_dists.tiff", p_1a, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.5, dpi = 600) # paper="letter"




# Show that absolute of normal samples generated under null hypothesis follow folded normal
#-------------------------------------------------------------------------

# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data                |
n_trials <- 1e3
mus <- rep(0, n_trials)#mrunif(1e4, min = -5, max = 5)
n_obs <- 50

# Generate random samples based on random mu values
x_norm <- t(sapply(mus, function (x) rnorm(n_obs, mean = x, sd = 1),simplify = TRUE))
x_abs_norm <- t(sapply(mus, function (x) abs(rnorm(n_obs, mean = x, sd = 1)),simplify = TRUE))
x_fnorm <- t(sapply(mus, function (x) rfoldnorm(n_obs, mean = x, sd = 1),simplify = TRUE))

# Test if norm is different from abs_norm
p_value_norm2abs_norm <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                                as.vector(x_abs_norm[i,]), alternative = "two.sided", exact = TRUE, 
                                tol=1e-8, simulate.p.value=FALSE)$p)

# Test if norm is different from folded norm
p_value_norm2fnorm <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_norm[i,]), 
                             as.vector(x_fnorm[i,]), alternative = "two.sided", exact = TRUE, 
                             tol=1e-8, simulate.p.value=FALSE)$p)

# Test if abs_norm is different from folded norm
p_value_abs_norm2fnorm <- sapply(1:nrow(x_norm), function(i) ks.test(as.vector(x_abs_norm[i,]),
                                 as.vector(x_fnorm[i,]), alternative = "two.sided", exact = TRUE, 
                                 tol=1e-8, simulate.p.value=FALSE)$p)
df_1b <- rbind(tibble(d = factor("N : |N|"), x = p_value_norm2abs_norm),
               tibble(d = factor("N : F[N]"), x = p_value_norm2fnorm), 
               tibble(d = factor("|N| : F[N]"), x = p_value_abs_norm2fnorm))

p_1b <- ggplot(df_1b, aes(x=d, y=x)) + 
  geom_violin(scale="width",bw=.05,kernel="gaussian") +
  xlab("Sample Comparisons") + ylab("K.S. P Value") + geom_boxplot(width=0.1,outlier.shape = NA) +
  theme_small_classic()
p_1b
# Export to TIF
save_plot("figure/figure_1b_dist_comparison.tiff", p_1b, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600) # paper="letter"




# Show % of trials that have significant p-values
# Under null disctribution it would be expected to equal alpha
bin_n_an <- prop.test(sum(p_value_norm2abs_norm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)
bin_n_fn <- prop.test(sum(p_value_norm2fnorm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)
bin_an_fn <- prop.test(sum(p_value_abs_norm2fnorm<0.05), n_trials, p=0.05, conf.level=1-0.05/3, correct = FALSE)

df_1c <- tibble(d = factor(c("N : |N|","N : F[N]","|N| : F[N]")), 
       estimate = c(bin_n_an$estimate,bin_n_fn$estimate,bin_an_fn$estimate),
       lcl = c(bin_n_an$conf.int[1],bin_n_fn$conf.int[1],bin_an_fn$conf.int[1]),
       ucl = c(bin_n_an$conf.int[2],bin_n_fn$conf.int[2],bin_an_fn$conf.int[2]),
       p_value = c(bin_n_an$p.value,bin_n_fn$p.value,bin_an_fn$p.value))

p_1c <- ggplot(df_1c, aes(x=d, y=estimate)) + 
  geom_hline(yintercept = 0.05, linetype="dashed", color="red", size=.5,alpha=.2) +
  geom_hline(yintercept = 1, linetype="dashed", color="black", size=.5,alpha=.2) +
  geom_point(size=.1) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2) +
  xlab("Sample Comparisons") +  ylab("Fraction Trials p<a  ") +
  theme_small_classic()
p_1c  
save_plot("figure/figure_1c_dist_comparison.tiff", p_1c, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 2, dpi = 600) # paper="letter"




# Plot difference between P(N) and P(F[N]) over mu
sort_x_abs_norm = sort(x_abs_norm[1,])
sort_x_fnorm = sort(x_fnorm[1,])

# Fit line to difference and show 95% CI is equal to zero
# Fit line (not used) to show correspondence
fit <- lm(formula = sort_x_abs_norm ~ 0 +sort_x_fnorm)
confint(fit,level=0.95)
# Define data for plotting
df_1d = tibble(x=sort_x_fnorm, y=sort_x_abs_norm)

p_1d <- ggplot(df_1d, aes(x=x, y=y)) + 
  geom_point(size=0.5) + xlab("Folded Normal Sample") + ylab("Absolute Normal Sample   ") +
  geom_abline(intercept = 0, slope = 1, size=0.25) +
  theme_small_classic()
p_1d
save_plot("figure/figure_1d_qq_plot.tiff", p_1d, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) # paper="letter"






# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
sd = 1
n_trials = 1e4
n_samples = 50
mus = seq(-3,3, by = 1)

# Simulation results to be calculated
# Rows: simulation #
# Cols: swept mu value
samples_mu <- matrix(rep(0,n_trials*length(mus)), ncol = length(mus))
samples_sd <- matrix(rep(0,n_trials*length(mus)), ncol = length(mus))
for (n in seq(1,length(mus),1)) {
  # Generate  samples for this simulation, calculate mean and sd for each, store in column
  obs_data <- abs(matrix(rnorm(n_trials*n_samples,mean = mus[n], sd=1), nrow = n_trials, byrow = TRUE))
  # Mean of samples for each trial
  samples_mu[,n] <- t(apply(obs_data,1,mean))
  # Standard deviation of samples for each trial
  samples_sd[,n] <- t(apply(obs_data,1,sd))
}

# Collapse matrix of mu's and sd's into a data frame
tbl_mu <-as.table(samples_mu,keep.rownames=FALSE,row.names=NULL,responseName=y)
colnames(tbl_mu) <- factor(mus)
df_mu <- as.data.frame(tbl_mu) %>% rename(index=Var1,mu = Var2, sample_mean = Freq)
tbl_sd <-as.table(samples_sd, keep.rownames = FALSE,row.names=NULL,responseName=y)
colnames(tbl_sd) <- factor(mus)
df_sd <- as.data.frame(tbl_sd) %>% rename(index=Var1,mu = Var2, sample_sd = Freq)

# Reformatting datatable 
df_results <- df_mu
df_results$sample_sd <- df_sd$sample_sd
df_results$mu <- as.factor(df_results$mu)
df_results$sample_sd_sig <- '*'
df_results$sample_mean_sig <-'*'

# Statistics of groups
# Compute the analysis of variance
res.aov <- aov(sample_mean ~ mu, data = df_results)
# Summary of the analysis
summary(res.aov)
pairwise.t.test(df_results$sample_mean, df_results$mu, p.adj = "bonf")
pairwise.t.test(df_results$sample_sd, df_results$mu, p.adj = "bonf")

# Sample mean versus population mean
p_1e <- ggplot(df_results, aes(x=mu, y=sample_mean)) + 
  geom_violin(color="blue",lwd = .2) + 
  #annotate("text",x=factor(mus),y=4,label='*')+
  geom_text(data=data.frame(), aes(x = factor(mus), y=rep(4,length(mus)), label=rep(c('*'), length(mus)) ), size=6) +
  coord_cartesian(clip = 'off') +
  geom_boxplot(width=.2, outlier.shape = NA,lwd = .2) +
  xlab("Population Mean") + ylab("Sample Mean") +
  theme_small_classic()
p_1e
save_plot("figure/figure_1e_mean_changes_sample_mean.tiff", p_1e, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600) # paper="letter"

# Sample sd versus population mean
p_1f <- ggplot(df_results, aes(x=mu, y=sample_sd)) + 
  geom_violin(color="blue", lwd = .2) + 
  geom_text(data=data.frame(), aes(x = factor(mus), y=rep(1.5,length(mus)), label=rep(c('*'), length(mus)) ), size=6) +
  coord_cartesian(clip = 'off') +
  geom_boxplot(width=.2, outlier.shape = NA, lwd = .2) +
  xlab("Population SD") + ylab("Sample SD") +
  theme_small_classic()
p_1f
save_plot("figure/figure_1f_mean_changes_sample_sd.tiff", p_1f, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600) # paper="letter"







# With fixed mean, change the sd from .1:.1:2, calculated ffolded mean and sd, and show they vary


# Show how mean is changed as mu increases from zero
# With fixed sd, sweep mean from -3:.25:3, calculate folded mean and sd, show they vary
pop_sds = seq(.5,5, by = .5)
n_trials = 1e4
n_samples = 50
mu = 0

# Simulation results to be calculated
# Rows: simulation #
# Cols: swept mu value
samples_mu <- matrix(rep(0,n_trials*length(pop_sds)), ncol = length(pop_sds))
samples_sd <- matrix(rep(0,n_trials*length(pop_sds)), ncol = length(pop_sds))
for (n in seq(1,length(pop_sds),1)) {
  # Generate  samples for this simulation, calculate mean and sd for each, store in column
  obs_data <- abs(matrix(rnorm(n_trials*n_samples,mean = mu, sd=pop_sds[n]), nrow = n_trials, byrow = TRUE))
  # Mean of samples for each trial
  samples_mu[,n] <- t(apply(obs_data,1,mean))
  # Standard deviation of samples for each trial
  samples_sd[,n] <- t(apply(obs_data,1,sd))
}

# Collapse matrix of mu's and sd's into a data frame
tbl_mu <-as.table(samples_mu,keep.rownames=FALSE,row.names=NULL,responseName=y)
colnames(tbl_mu) <- factor(pop_sds)
df_mu <- as.data.frame(tbl_mu) %>% rename(index=Var1,pop_sd = Var2, sample_mean = Freq)
tbl_sd <-as.table(samples_sd, keep.rownames = FALSE,row.names=NULL,responseName=y)
colnames(tbl_sd) <- factor(sds)
df_sd <- as.data.frame(tbl_sd) %>% rename(index=Var1,pop_sd = Var2, sample_sd = Freq)

# Reformatting datatable 
df_results <- df_mu
df_results$sample_sd <- df_sd$sample_sd
df_results$pop_sd <- as.factor(df_results$pop_sd)
df_results$sample_sd_sig <- '*'
df_results$sample_mean_sig <-'*'

# Statistics of groups
# Compute the analysis of variance
summary(aov(sample_mean ~ pop_sd, data = df_results))
summary(aov(sample_sd ~ pop_sd, data = df_results))
pairwise.t.test(df_results$sample_mean, df_results$pop_sd, p.adj = "bonf")
pairwise.t.test(df_results$sample_sd, df_results$pop_sd, p.adj = "bonf")


plot_label_height <- function (var,line_num, pad_mult) {
  line_ypos <- pad_mult*(max(var)-min(var)) + max(var) + line_num * pad_mult/2*(max(var)-min(var))
}

# Sample mean versus population mean
p_1g <- ggplot(df_results, aes(x=pop_sd, y=sample_mean)) + 
  geom_violin(color="blue",lwd = .2) + 
  #annotate("text",x=factor(mus),y=4,label='*')+
  geom_text(data=data.frame(), aes(x = factor(pop_sds), y=rep(3,length(pop_sds)), label=sample_mean_sig ), size=6) +
  coord_cartesian(clip = 'off') +
  geom_boxplot(width=.2, outlier.shape = NA,lwd = .2) +
  xlab("Population SD") + ylab("Sample Mean") +
  theme_small_classic()
p_1g
save_plot("figure/figure_1e_mean_changes_sample_mean.tiff", p_1g, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600) # paper="letter"

# Sample sd versus population mean
p_1f <- ggplot(df_results, aes(x=pop_sd, y=sample_sd)) + 
  geom_violin(color="blue", lwd = .2) + 
  geom_text(data=data.frame(), aes(x = factor(pop_sds), y=rep(plot_label_height(sample_sd),length(pop_sds)), label=sample_sd_sig ), size=6) +
  coord_cartesian(clip = 'off') +
  geom_boxplot(width=.2, outlier.shape = NA, lwd = .2) +
  xlab("Population SD") + ylab("Sample SD") +
  theme_small_classic()
p_1f
save_plot("figure/figure_1f_mean_changes_sample_sd.tiff", p_1f, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3, dpi = 600) # paper="letter"

