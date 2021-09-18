 

#' A series of examples to illustrate the three parameters of unscaled agreement
#' and relative agreement. Each example compares the difference in means or 
#' relative difference in means between a hypothetical control and experiment
#' group.


# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(ggplot2)
p_load(dplyr)
p_load(cowplot)
p_load(RColorBrewer)
p_load(stringr)


# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
fig_num = "0" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(file.path(getwd(), fig_path), showWarnings = FALSE)
grp_color = brewer.pal(3, "Set1")
yaxis_font_size = 7
base_font_size = 8

# Figrue sizes, height x length
ggsize = c(1.2, .65)

# User defined functions
plot_experiments <- function(xa1,xb1,xa2,xb2, fig_num, fig_path, base_name, 
                             data_scale = "absolute", ylims = c(NULL,NULL), 
                             dm_ylims = c(NULL,NULL), dm_ybreaks = c(NULL,NULL), 
                             conf.level.1 = 0.95, conf.level.2 = 0.95) {
  save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  #' @description 
  #' 
  #' 
  #' @param xa1 vector of measurements from group a, experiment 1
  #' @param xb1 vector of measurements from group a, experiment 1
  #' @param xa2 vector of measurements from group a, experiment 1
  #' @param xb2 vector of measurements from group a, experiment 1
  #' @param fig_num string of number label used for basename of all figures
  #' @param fig_path path to export figures to disk
  #' @param base_name base name for exported figures
  #' @param data_scale scale either "absolute" or "relative"
  #' @param ylims limits of y axis scale for raw data
  #' @param dm_ylims limits of y axis of difference in means
  #' @param dm_ybreaks vector of breakpoints for y axis of difference in means
  #' 
  #' @return null exports various plots to disk

  # load(file = "temp/debug.RData")
  
  na1 = length(xa1); nb1 = length(xb1)
  na2 = length(xa2); nb2 = length(xb2)  
  blank_df <- tibble(group=c("a","a"), y=ylims)
  
  # Experiment 1
  df <- tibble(group = as.factor(c(rep("a",length(xa1)), rep("b",length(xb1)))), y = c(xa1,xb1))
  df_stats <- df %>% group_by(group) %>% summarise(mean = mean(y),sd = sd(y))
  gg1 <- ggplot(data=df, aes(x=group, y=y)) +
    geom_jitter(width = 0.2, size=0.6, alpha = 0.5,color = grp_color[1],shape=16) + 
    geom_point(data = df_stats, aes(x=group, y=mean), size=6, shape='-') +
    geom_blank(data=blank_df, aes(x=group, y=y)) + coord_cartesian(ylim = ylims) +
    theme_classic(base_size=yaxis_font_size) + xlab("Exp 1") + ylab("Data") +
    geom_errorbar(data = df_stats, aes(x=group,ymin=mean-sd, ymax=mean+sd, y=mean), 
                  width=0.3) +
    theme(axis.text=element_text(size=base_font_size),
          axis.text.y = element_text(size=yaxis_font_size),
          axis.title.x = element_text(color = grp_color[1]))
  gg1
  # if (max(df$y)<0.02) {gg1 <- gg1 + scale_y_continuous(labels = scales::scientific)}
  save_plot(paste(fig_path, "/F",fig_num, "_", base_name,"_E1", ".tiff", sep = ""),
            gg1, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])
  
  
  # Experiment 2
  df <- tibble(group = as.factor(c(rep("a",length(xa2)),rep("b",length(xb2)))),y = c(xa2,xb2))
  df_stats <- df %>% group_by(group) %>% summarise(mean = mean(y),sd = sd(y))
  gg2 <- ggplot(data=df, aes(x=group, y=y)) +
    geom_jitter(width = 0.2, size=0.6, alpha = 0.5,color = grp_color[2],shape=16) +
    geom_point(data = df_stats, aes(x=group, y=mean), size=6, shape='-') +
    geom_blank(data=blank_df, aes(x=group, y=y)) + coord_cartesian(ylim = ylims) +
    theme_classic(base_size=yaxis_font_size) + xlab("Exp 2") + ylab("Data") +
    geom_errorbar(data = df_stats, aes(x=group,ymin=mean-sd, ymax=mean+sd, y=mean), 
                  width=0.4) +
    expand_limits(y = min(df$y) + 1.2*(max(df$y)-min(df$y))) +
    theme(axis.text=element_text(size=base_font_size),
          axis.text.y = element_text(size=yaxis_font_size),
          axis.title.x = element_text(color = grp_color[2]))
  gg2
  if (max(df$y)<0.02) {gg2 <- gg2 + scale_y_continuous(labels = scales::scientific)}
  save_plot(paste(fig_path, "/F", fig_num, "_", base_name,"_E2", ".tiff", sep = ""),
            gg2, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])
  
  # Difference in Means
  blank_df <- tibble(group=c(1,2), y=dm_ylims)
  # Stats for each group
  dfo <- tibble(mu_d1 = mean(xb1) - mean(xa1), sigma_d1 = sqrt(sd(xa1)^2 + sd(xb1)^2),
                mu_dm1 = mean(xb1) - mean(xa1), sigma_dm1 = sqrt(sd(xa1)^2/na1 + sd(xb1)^2/nb1),
                mu_d2 = mean(xb2) - mean(xa2), sigma_d2 = sqrt(sd(xa2)^2 + sd(xb2)^2),
                mu_dm2 = mean(xb2) - mean(xa2), sigma_dm2 = sqrt(sd(xa2)^2/na2 + sd(xb2)^2/nb2),
                rsigma_d1 = sigma_d1/(mean(xa1)+0.5*mu_d1),
                rsigma_d2 = sigma_d2/(mean(xa2)+0.5*mu_d2),
                rmu_dm1 = mu_dm1/mean(xa1), rsigma_dm1 = sigma_dm1/(mean(xa1) + 0.5*mu_d1),
                rmu_dm2 =    mu_dm2/mean(xa2), rsigma_dm2 = sigma_dm2/(mean(xa2) + 0.5*mu_d2),
                df1=length(xa1)+length(xb1)-2, df2=length(xa2)+length(xb2)-2)

  # n in this case is number of samples for distribution histogram
  n1 = 1e6; n2 = 1e6
  
  # Plot 95% CI of difference in means for both experiments
  if (data_scale=="relative") {
    # browser();
    
    grp1_ci95 = quantile(rnorm(1e5,mean = dfo$mu_dm1, sd = dfo$sigma_dm1)/mean(xa1), 
                         c((1-conf.level.1)/2, 1 - (1-conf.level.1)/2)) 
    grp2_ci95 = quantile(rnorm(1e5,mean = dfo$mu_dm2, sd = dfo$sigma_dm2)/mean(xa2), 
                         c((1-conf.level.2)/2, 1 - (1-conf.level.2)/2)) 
    
    df_ci <- tibble(group = factor(c(1,2)), y = c(dfo$mu_dm1/mean(xa1), dfo$mu_dm2/mean(xa2)), 
                    ymin = c(grp1_ci95[1], grp2_ci95[1]), 
                    ymax = c(grp1_ci95[2], grp2_ci95[2]))
    
    # Both experiments
    gg3 <- ggplot(data = df_ci, aes(x=group, y=y)) +
      geom_pointrange(aes(color=group, ymin = ymin, ymax = ymax),size=0.5,
                      fatten = 0.2, shape = 16) +
      xlab("Exp") + ylab("Relative DM") +
      coord_cartesian(ylim = dm_ylims) + #geom_blank(data=blank_df, aes(x=group, y=y)) +
      scale_color_manual(values=grp_color[1:2]) +
      theme_classic(base_size=7) +
      theme(axis.text=element_text(size=base_font_size),legend.position = "none",
            axis.text.y = element_text(size=yaxis_font_size),
            axis.text.x = element_text(colour = grp_color))
    gg3
    
  } else {y_label = expression(mu[DM]); denominator = 1; 
  
  # browser()
  
  
  grp1_ci95 = quantile(rnorm(1e5,mean = dfo$mu_dm1, sd = dfo$sigma_dm1), 
                       c((1-conf.level.1)/2, 1 - (1-conf.level.1)/2)) 
  grp2_ci95 = quantile(rnorm(1e5,mean = dfo$mu_dm2, sd = dfo$sigma_dm2), 
                       c((1-conf.level.2)/2, 1 - (1-conf.level.2)/2)) 

  df_ci <- tibble(group = factor(c(1,2)), y = c(dfo$mu_dm1, dfo$mu_dm2), 
                  ymin = c(grp1_ci95[1], grp2_ci95[1]), 
                  ymax = c(grp1_ci95[2], grp2_ci95[2]))
  
  # Both experiments
  gg3 <- ggplot(data = df_ci, aes(x=group, y=y)) +
    geom_pointrange(aes(color=group, ymin = ymin, ymax = ymax), size=0.5,
                    fatten = 0.2, shape = 16) +
    xlab("Exp") + ylab("DM") +
    geom_blank(data=blank_df, aes(x=group, y=y)) +
    scale_y_continuous(breaks = dm_ybreaks) +
    scale_color_manual(values=grp_color[1:2]) +
    theme_classic(base_size = 7) +
    theme(axis.text=element_text(size=base_font_size),
          legend.position = "none",
          axis.text.y = element_text(size=yaxis_font_size),
          axis.text.x = element_text(colour = grp_color))
  gg3
  }
  
  save_plot(paste(fig_path, "/F", fig_num, "_", base_name,"_D", ".tiff", sep = ""),
            gg3, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])
  print(t(dfo))
  
  ggplots = list(gg1, gg2, gg3);
  names(ggplots) <- c("gg1","gg2","gg3")
  return(ggplots)
  # browser();
}



#-------------------------------------------------------------------------------

# Absolute Scale

#-------------------------------------------------------------------------------

# Decreased mean offset for experiment 2
set.seed(0)
na1=10; nb1=10; na2=10; nb2=10;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)



# Same P, weaker practical significance
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(.5) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(.5) + 10
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(4^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(4^2/2) + 10
# Plot experiments and Difference of Means
gg_mu <- plot_experiments(xa1, xb1, xa2, xb2, fig_num=fig_num, fig_path=fig_path, base_name="unsc_mu", 
                 ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2))
sprintf("p-val) 1: %.2e, 2: %.2e", t.test(xa1, xb1)$p.value,t.test(xa2, xb2)$p.value)
sprintf("deltaM) 1: %.3f, 2: %.3f", rmdm_tdist(xa1,xb1),rmdm_tdist(xa2,xb2))




# larger p, weaker practical significance
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(.3) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(.3) + 10.4
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(3.6^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(3.6^2/2) + 10.4
# Plot experiments and Difference of Means
gg_mu <- plot_experiments(xa1, xb1, xa2, xb2, fig_num=fig_num, fig_path=fig_path, base_name="pval_larger", 
                          ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2))
sprintf("p-val) 1: %.2e, 2: %.2e", t.test(xa1, xb1)$p.value,t.test(xa2, xb2)$p.value)
sprintf("deltaM) 1: %.3f, 2: %.3f", rmdm_tdist(xa1,xb1),rmdm_tdist(xa2,xb2))



# Smaller p, weaker practical significance
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(3) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(3) + 10.1
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(3) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(3) + 11.5
# Plot experiments and Difference of Means
gg_mu <- plot_experiments(xa1, xb1, xa2, xb2, fig_num=fig_num, fig_path=fig_path, base_name="pval_smaller", 
                          ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2))
sprintf("p-val) 1: %.2e, 2: %.2e", t.test(xa1, xb1)$p.value,t.test(xa2, xb2)$p.value)
sprintf("deltaM) 1: %.3f, 2: %.3f", rmdm_tdist(xa1,xb1),rmdm_tdist(xa2,xb2))




