

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
title_font_size=9;
equ_font_size=8;

# Initial parameters
nsamples = 1e6
ngroups = 2
set.seed(0)
dist_text_size=3;
title_font_size=9;
equ_font_size=8;



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
  facet_grid(. ~ mu, scales = "free_x", labeller=label_parsed) +
  theme_classic() + ylab("Probability") +
  scale_fill_discrete(name = "", labels = c(expression(f(x)), expression(abs(phantom(0)*f(x)~phantom())))) + 
  theme(plot.title = element_text(hjust = 0.5, size = title_font_size,family="serif"),legend.position = "none",
        text = element_text(size=9))
# p_1a
# Export to TIF
save_plot("figure/figure_1a_example_dists.tiff", p_1a, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 3.5, dpi = 600) # paper="letter"




# Test to show that absolute of normal samples generated under null hypothesis follow folded normal
#-------------------------------------------------------------------------

# Randomly generate normal data for     
#     Normal data,                      |
#     Normal data, then apply ABS()     |   x   1000
#     Folded normal data                |
mus <- rep(0, 1e4)#mrunif(1e4, min = -5, max = 5)
n_obs <- 200

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


df_1b <- rbind(tibble(d = "N : |N|", x = p_value_norm2abs_norm),
               tibble(d = "N : F[N]", x = p_value_norm2fnorm), 
               tibble(d = "|N| : F[N]", x = p_value_abs_norm2fnorm))
 
p_1b <- ggplot(df_1b, aes(x=d, y=x)) + 
  geom_violin(scale="width",bw=.05,kernel="gaussian") +
  xlab("Comparison") + ylab("K.S. P Value") + geom_boxplot(width=0.1,outlier.shape = NA)
p_1b


a <- hist(p_value_abs_norm2fnorm, breaks=9)

chisq.test(hist(p_value_abs_norm2fnorm))

library(spgs)
chisq.unif.test(hist(p_value_abs_norm2fnorm), interval=c(0,1))
                
                
uniform.test(hist(p_value_abs_norm2fnorm), B = 1000)

ks.test(p_value_abs_norm2fnorm,"punif",0,1,exact=TRUE, simulate.p.value=TRUE)

# Plot difference between P(N) and P(F[N]) over mu


# Fit line to difference and show 95% CI is equal to zero




# Quantile plot that shows samples follows folded normal



# Show how mean is changed as x_bar increases from zero



# Show how std increases as x_bar changes from zero