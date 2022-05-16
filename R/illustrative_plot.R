

if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(ggplot2)
p_load(dplyr)
p_load(cowplot)
p_load(RColorBrewer)
p_load(stringr)
source("R/aces.R")


# User defined functions
plot_experiment_pair <- function(xa1,xb1,xa2,xb2, fig_num, fig_path, base_name, 
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
  blank_df <- tibble(group=c("X","Y"), y=ylims)
  
  # Experiment 1
  df <- tibble(group = as.factor(c(rep("X",length(xa1)), rep("Y",length(xb1)))), y = c(xa1,xb1))
  df_stats <- df %>% group_by(group) %>% summarise(mean = mean(y),sd = sd(y))
  gg1 <- ggplot(data=df, aes(x=group, y=y)) +
    geom_jitter(width = 0.2, size=0.6, alpha = 0.5,color = grp_color[1],shape=16) + 
    geom_point(data = df_stats, aes(x=group, y=mean), size=7, shape='-') +
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
  df <- tibble(group = as.factor(c(rep("X",length(xa2)),rep("Y",length(xb2)))),y = c(xa2,xb2))
  df_stats <- df %>% group_by(group) %>% summarise(mean = mean(y),sd = sd(y))
  gg2 <- ggplot(data=df, aes(x=group, y=y)) +
    geom_jitter(width = 0.2, size=0.6, alpha = 0.5,color = grp_color[2],shape=16) +
    geom_point(data = df_stats, aes(x=group, y=mean), size=7, shape='-') +
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
    
    rci1 <- credint(xa1,xb1,conf.level = conf.level.1, rand.seed = 0, relative = TRUE)
    rci2 <- credint(xa2,xb2,conf.level = conf.level.2, rand.seed = 0, relative = TRUE)
    df_ci <- tibble(group = factor(c(1,2)), y = c(dfo$mu_dm1/mean(xa1), dfo$mu_dm2/mean(xa2)), 
                    ymin = c(rci1[1], rci2[1]), 
                    ymax = c(rci1[2], rci2[2]))
    
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
  
  
  rci1 <- credint(xa1,xb1,conf.level = conf.level.1, rand.seed = 0, relative = FALSE)
  rci2 <- credint(xa2,xb2,conf.level = conf.level.2, rand.seed = 0, relative = FALSE)
  df_ci <- tibble(group = factor(c(1,2)), y = c(dfo$mu_dm1, dfo$mu_dm2), 
                  ymin = c(rci1[1], rci2[1]), 
                  ymax = c(rci1[2], rci2[2]))
  
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
