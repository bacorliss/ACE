 



# p.adjust(p, method = p.adjust.methods, n = length(p))


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
library(binom)
library(VGAM)
library(gplots)
library(RColorBrewer)
# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
library(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
library(TOSTER)



# choose colors for plotting
color_pal = brewer.pal(4, "Set1")
fig_basename="f_0"

# Parameters

nsamples = 1e4


# Generate original data
# Data used for all example figues
rand_seed <- 0
xa <- rnorm(n = 10, mean = 10, sd = 1)
xb <- rnorm(n = 10, mean = 10, sd = 1)

# Normalize SD
xa <- (xa-10)/sd(xa)+10
xb <- (xb-10)/sd(xb)+10
# Normalize mean
xa <- xa - mean(xa)+10
xb <- xb - mean(xb)+10

xa_bar <- mean(xa)
xa_sd  <- sd(xa)
xb_bar <- mean(xb)
xb_sd   <- sd(xb)
rank_str <- "test"
is_exported <- TRUE



plot_null_result <- function(xa, xa_bar, xa_sd, xb, xb_bar,xb_sd, rank_str,is_exported) {
  
  # Transform data
  txa <- xa_sd * (xa - mean(xa)) + +mean(xa) + xa_bar
  txb <- xb_sd * (xb - mean(xb)) + +mean(xb) + xb_bar
  # Assemble into data frame
  df1 <- rbind(tibble(grp = "Control", y = txa), tibble(grp = "Tx", y = txb))
  
  # Get standard deviation and mean for error bars
  sdf1 <- df1 %>% group_by(grp) %>% summarize(mean = mean(y), sd = sd(y))
  
  # Plotting
  p_1 <- ggplot(data=df1, aes(x = as.factor(grp), y = y)) +
    geom_point(position = position_jitter(w = 0.15, h = 0), color="grey", size=0.5) +
    geom_point(data = sdf1, aes(x = grp, y = mean), size=1) +
    geom_errorbar(data = sdf1, aes(x = grp, y = mean, ymin = mean - sd, ymax = mean + sd), width = 0.5) + 
    theme_classic(base_size=8) + theme(legend.position="none")+
    coord_cartesian(ylim=c(0,20)) + xlab("") + ylab("Metric")
  print(p_1)
  # Stats
  tt <- t.test(subset(df1, grp=="Tx")$y, subset(df1, grp=="Control")$y)
  
  # Print stats
  print(do.call(sprintf, c('Mean: [%.1f, %.1f]', as.list(c(mean(txa),mean(txb))))))
  print(do.call(sprintf, c('STD: [%.1f, %.1f]', as.list(c(sd(txa),sd(txb))))))
  print(do.call(sprintf, c('DOM: %.1f', as.list(mean(txb)-mean(txa)))))
  print(do.call(sprintf, c('SDpooled: %.1f', as.list(sqrt(sd(txa)^2+sd(txb)^2)))))
  print(do.call(sprintf, c('CI: [%.1f, %.1f]', as.list(unname(round(tt$conf.int,2)+10)))))
  print(do.call(sprintf, c('R-CI: [%.1f, %.1f]', as.list(unname(round(tt$conf.int/10,3)*100)))))
  print(do.call(sprintf, c('pval: %.3f', as.list(unname(tt$p.value)))))
  
  # Export
  if (is_exported) {
  save_plot( paste("figure/", 
          fig_basename, sprintf("_%s.tiff", rank_str), 
          sep = ""),
          p_1, ncol = 1, nrow = 1, base_asp = 3, dpi = 600,  
          base_height = 1.5, base_width = 1.5)
  }
  
}

# Easy Example
plot_null_result(xa, 0, 1, xb, 0.2, 1, rank_str = "r1_1",is_exported = TRUE)
plot_null_result(xa,0, 4,xb, 0.8, 4,rank_str = "r1_2",is_exported = TRUE)



# Hard Example
plot_null_result(xa, 0, 5,xb, 1.5, 5, rank_str = "r2_1", is_exported = TRUE)

plot_null_result(xa, +1, 3,xb, -1, 3, rank_str = "r2_2", is_exported = TRUE)

plot_null_result(xa, -1, 2,xb, +3, 7, rank_str = "r2_3", is_exported = TRUE)

plot_null_result(xa, 0, 7,xb, 0.1, 7, rank_str = "r2_4", is_exported = TRUE)

plot_null_result(xa, 0, 1.5,xb, 1, 1.5, rank_str = "r2_5", is_exported = TRUE)

