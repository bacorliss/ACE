 


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
p_load(binom)
p_load(VGAM)
p_load(gplots)
p_load(RColorBrewer)
# https://cran.r-project.org/web/packages/equivalence/equivalence.pdf
p_load(equivalence)
# https://cran.rstudio.com/web/packages/TOSTER/vignettes/IntroductionToTOSTER.html
p_load(TOSTER)


source("R/mdm.r")


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

xa_shift <- mean(xa)
xa_spread  <- sd(xa)
xb_shift <- mean(xb)
xb_spread   <- sd(xb)
rank_str <- "test"
is_exported <- TRUE
y_label="Metric"



plot_null_result <- function(xa, xa_shift, xa_spread, xb, xb_shift,xb_spread, rank_str,is_exported, y_label="Metric") {
  # [xa]:  data [xa_shift]: how much top shift mean [xa_spread]: scale factor for sd
  # [xb]:  data [xb_shift]: how much top shift mean [xb_spread]: scale factor for sd
  
  # Transform data
  txa <- xa_spread * (xa - mean(xa)) + +mean(xa) + xa_shift
  txb <- xb_spread * (xb - mean(xb)) + +mean(xb) + xb_shift
  # Assemble into data frame
  df1 <- rbind(tibble(grp = "Control", y = txa), tibble(grp = "Tx", y = txb))
  
  # Get standard deviation and mean for error bars
  sdf1 <- df1 %>% group_by(grp) %>% summarize(mean = mean(y), sd = sd(y))
  
  # Plotting
  p_1 <- ggplot(data=df1, aes(x = as.factor(grp), y = y)) +
    geom_point(position = position_jitter(w = 0.15, h = 0), color="grey", size=0.5) +
    geom_point(data = sdf1, aes(x = grp, y = mean), size=1) +
    geom_errorbar(data = sdf1, aes(x = grp, y = mean, ymin = mean - sd, ymax = mean + sd), width = 0.5) + 
    theme_classic(base_size=8) + theme(legend.position="none") +
    xlab("") + ylab(y_label) #+ coord_cartesian(ylim=c(0,20))
  print(p_1)
  # Stats
  tt <- t.test(subset(df1, grp=="Tx")$y, subset(df1, grp=="Control")$y, 
               paired = FALSE, var.equal = FALSE,conf.level = 0.95)
  print(tt)
  mhd <- mhd_2sample_paired(txa, txb, alpha = 0.05)
  
  # Print stats
  print(do.call(sprintf, c('x_bar: [%.1f, %.1f]', as.list(c(mean(txa),mean(txb))))))
  print(do.call(sprintf, c('s: [%.1f, %.1f]', as.list(c(sd(txa),sd(txb))))))
  print(do.call(sprintf, c('pval: %.3f', as.list(unname(tt$p.value)))))
  print(do.call(sprintf, c('x_bar[D]: %.1f', as.list(mean(txb)-mean(txa)))))
  #print(do.call(sprintf, c('x_bar[RD]: %.1f', as.list((mean(txb)-mean(txa))/mean(txa)))))
  print(do.call(sprintf, c('s[D]: %.1f', as.list(sqrt(sd(txa)^2+sd(txb)^2)))))
  #print(do.call(sprintf, c('cv[D]: %.1f', as.list(sqrt(sd(txa)^2+sd(txb)^2)/mean(txa)))))
  
  
  print(do.call(sprintf, c('CI: [%.1f, %.1f]', as.list(unname(round(tt$conf.int,2))))))
  print(do.call(sprintf, c('R-CI: [%.1f, %.1f]', as.list(unname(round(tt$conf.int/tt$estimate[2],3)*100)))))

  print(do.call(sprintf, c('MHD: %.1f', as.list(mhd))))
  print(do.call(sprintf, c('RMHD: %.1f', as.list(100*mhd/mean(txa)))))
  
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
plot_null_result(xa, 0, 5,xb,        1.5, 5, rank_str = "r2_1", is_exported = TRUE, y_label= "Metric 1")

plot_null_result(xa, 1*10+1, 1.5*3,xb, 10-1, 1.5*3, rank_str = "r2_2", is_exported = TRUE, y_label= "Metric 2")

plot_null_result(xa, 5*10-2, 5*2,xb, 5*10+3, 5*6, rank_str = "r2_3", is_exported = TRUE, y_label= "Metric 3")

plot_null_result(xa, 2+0, 5,xb,      2+0.1, 5, rank_str = "r2_4", is_exported = TRUE, y_label= "Metric 4")

plot_null_result(xa, +5, 3,xb,     7, 4, rank_str = "r2_5", is_exported = TRUE, y_label= "Metric 5")



# Magnitude of effecet size example

plot_null_result(xa,0, 4,xb, 2, 4,rank_str = "r3_1",is_exported = TRUE)
plot_null_result(xa,0, 4,xb, -2, 4,rank_str = "r3_2",is_exported = TRUE)

