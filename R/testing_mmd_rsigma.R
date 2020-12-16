

# Load package manager
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
# Load packages
p_load(ggplot2)
p_load(tibble)
p_load(RColorBrewer)
p_load(broom)
# p_load(gridExtra)
# p_load(grid)
p_load(tidyr)
p_load(cowplot)
p_load(dplyr)
# p_load(effsize)
p_load(boot)
# User defined libraries
source("mdm.R")
source("effect_size_contests.R")



n_obs = 50;


# Generate two samples corrected to exactly standard normal samples
#  (Mean shift to zero, sd scaled to 1)
# 
set.seed(0)
x_a_init = rnorm(n_obs, mean = 0, sd = 1)
x_b_init = rnorm(n_obs, mean = 0, sd = 1)
x_d_init = rnorm(n_obs, mean = 0, sd = 1)

# Normalized samples (post generation)
x_a = (x_a_init-mean(x_a_init))/sd(x_a_init)
x_b = (x_b_init-mean(x_b_init))/sd(x_b_init)
x_d = (x_d_init-mean(x_d_init))/sd(x_d_init)


Ca = 20;

u_1a = Ca*10

u_1d = 10

sd_d <- sapply(Ca * seq(.1, 1, by=0.1), function(s) sd(s * (x_d) + u_1d))
sd_md = sd_d/sqrt(n_obs)

mdm_d <- sapply(Ca * seq(.1, 1, by=0.1), function(s) mdm_normal(s * (x_d) + u_1d,
                                                           paired = FALSE))

mdm_d/u_1a
1.64*sd_md/u_1a



sd_md
mdm_d

(mdm_d-u_1d)/(1.645*(sd_md))



plot(2*sd_md, mdm_d)



rsd_x <- sd_x/1
rmdm_x <- mdm_x/1


plot(rsd_x,rmdm_x)



u_1a = 10
sigma_1a = 1
u_1d = 50
sigma_1d = 50
rsigma_1md = (sigma_1d/sqrt(n_obs))/u_1a



u_2a = 30
sigma_2a = 1
u_2d = 50
sigma_2d = 150;
rsigma_2md = (sigma_2d/sqrt(n_obs))/u_2a





x_1at <- x_a * sigma_1a + u_1a
x_1bt <- x_b * sqrt(sigma_1a^2 + sigma_1d^2) + (u_1a + u_1d)

x_2at <- x_a * sigma_2a + u_2a
x_2bt <- x_b * sqrt(sigma_2a^2 + sigma_2d^2) + (u_2a + u_2d)


mdm1 = mdm_normal(x_1at, x_1bt,paired = FALSE)
rmdm1 = mdm1/mean(x_1at)


mdm2 = mdm_normal(x_2at, x_2bt,paired = FALSE)
rmdm2 = mdm2/mean(x_2at)

