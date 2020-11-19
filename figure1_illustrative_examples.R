 


library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(stringr)


grp_color = brewer.pal(3, "Set1")

yaxis_font_size = 7
base_font_size = 8

# Figrue sizes, height x length
ggsize = c(1.2, .75)

# User defined functions
plot_experiments <- function(xa1,xb1,xa2,xb2, fig_num, fig_path, base_name, 
                             data_scale = "absolute", ylims = c(NULL,NULL), 
                             dm_ylims = c(NULL,NULL), dm_ybreaks = c(NULL,NULL)) {
  save(list = ls(all.names = TRUE), file = "temp/debug.RData",envir = environment())
  # load(file = "temp/debug.RData")
  
  na1 = length(xa1); nb1 = length(xb1)
  na2 = length(xa2); nb2 = length(xb2)  
  blank_df <- tibble(group=c("A","A"), y=ylims)
  
  # Experiment 1
  df <- tibble(group = as.factor(c(rep("A",length(xa1)), rep("B",length(xb1)))), y = c(xa1,xb1))
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
  save_plot(paste(fig_path, "F",fig_num, "_", base_name,"_E1", ".tiff", sep = ""),
            gg1, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])
  # Experiment 2
  df <- tibble(group = as.factor(c(rep("A",length(xa2)),rep("B",length(xb2)))),y = c(xa2,xb2))
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
  # if (max(c(xa1,xb1))<1e-2) {
    # gg2 <- gg2 + scale_y_continuous(labels=function(x){
      # format(str_replace(fshape=21ormatC(x, format = "e", digits = 1), 'e-0','e-'))})
  # }
  gg2
  if (max(df$y)<0.02) {gg2 <- gg2 + scale_y_continuous(labels = scales::scientific)}
  save_plot(paste(fig_path, "F", fig_num, "_", base_name,"_E2", ".tiff", sep = ""),
            gg2, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])
  
  # Difference in Means
  blank_df <- tibble(group=c(1,1), y=dm_ylims)
  # Stats for each group
  dfo <- tibble(mu_d1 = mean(xb1) - mean(xa1),sigma_d1 = sqrt(sd(xa1)^2 + sd(xb1)^2),
                mu_dm1 = mean(xb1) - mean(xa1), sigma_dm1 = sqrt(sd(xa1)^2/na1 + sd(xb1)^2/nb1),
                mu_d2 = mean(xb2) - mean(xa2), sigma_d2 = sqrt(sd(xa2)^2 + sd(xb2)^2),
                mu_dm2 = mean(xb2) - mean(xa2), sigma_dm2 = sqrt(sd(xa2)^2/na2 + sd(xb2)^2/nb2),
                rsigma_d1 = sigma_d1/(mean(xa1)+0.5*mu_d1),
                rsigma_d2 = sigma_d2/(mean(xa2)+0.5*mu_d2),
                rmu_dm1 = mu_dm1/mean(xa1), rsigma_dm1 = sigma_dm1/(mean(xa1) + 0.5*mu_d1),
                rmu_dm2 =    mu_dm2/mean(xa2), rsigma_dm2 = sigma_dm2/(mean(xa2) + 0.5*mu_d2),
                df1=length(xa1)+length(xb1)-2, df2=length(xa2)+length(xb2)-2)

  # n in this case is number of samples for distribution histogram
  n1 = 1e5; n2 = 1e5
  if (data_scale=="relative") {
    df_sample <- tibble(group = as.factor(c(rep(1,n1),rep(2,n2))), 
                        y = c(rnorm(n1,mean = dfo$mu_dm1, sd = dfo$sigma_dm1)/mean(xa1),
                              rnorm(n2,mean = dfo$mu_dm2, sd = dfo$sigma_dm2)/mean(xa2)))
    # Both experiments
    gg3 <- ggplot(data = df_sample, aes(x=group, y=y)) +
      geom_violin(aes(fill=group, width = 0.8), scale = "width", alpha=0.7) + theme_classic(base_size=7) + 
      xlab("Exp") + ylab("Rel. Diff. in Means    ") +
      geom_blank(data=blank_df, aes(x=group, y=y)) + coord_cartesian(ylim = dm_ylims) +
      scale_fill_manual(values=grp_color[1:2]) +
      theme(axis.text=element_text(size=base_font_size),legend.position = "none",
            axis.text.y = element_text(size=yaxis_font_size),
            axis.text.x = element_text(colour = grp_color))
    gg3
    
  } else {y_label = expression(mu[DM]); denominator = 1; 
  df_sample <- tibble(group = as.factor(c(rep(1,n1),rep(2,n2))), 
                      y = c(rnorm(n1,mean = dfo$mu_dm1, sd = dfo$sigma_dm1)/1,
                            rnorm(n2,mean = dfo$mu_dm2, sd = dfo$sigma_dm2)/1))
  # Both experiments
  gg3 <- ggplot(data = df_sample, aes(x=group, y=y)) +
    geom_violin(aes(fill=group,width = 0.8), scale = "width", alpha=0.7) + theme_classic(base_size=7) + 
    xlab("Exp") + ylab("Diff. in Means") +
    geom_blank(data=blank_df, aes(x=group, y=y)) +
    scale_y_continuous(breaks = dm_ybreaks) +
    scale_fill_manual(values=grp_color[1:2]) +
    theme(axis.text=element_text(size=base_font_size),legend.position = "none",
          axis.text.y = element_text(size=yaxis_font_size),
          axis.text.x = element_text(colour = grp_color))
  gg3
  }
  
  save_plot(paste(fig_path, "F", fig_num, "_", base_name,"_D", ".tiff", sep = ""),
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

fig_num = "1" 
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/F",fig_num, "/",sep="")

# Decreased mean offset for experiment 2
rand_seed <- 0
na1=5; nb1=5; na2=5; nb2=5;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)

# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(0.9^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(0.9^2/2) + 10.5
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(0.9^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(0.9^2/2) + 12
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2,fig_num=fig_num, fig_path=fig_path, base_name="mu_box", 
                 ylims = c(6,14), dm_ylims = c(-3.7,4), dm_ybreaks=seq(-3,3,3))


# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1.9^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1.9^2/2) + 10.5
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, base_name="sigma_box", 
                 ylims = c(6,14), dm_ylims = c(-3.7,4), dm_ybreaks=seq(-3,3,3))


# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na2=35; nb2=35;
a3 <- rnorm(na2, mean = 0, sd = 1)
b3 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa3 <- (a3-mean(a3))/sd(a3)*sqrt(1^2/2) + 10
xb3 <- (b3-mean(b3))/sd(b3)*sqrt(1^2/2) + 10.5
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa3,xb3, fig_num=fig_num, fig_path=fig_path, base_name="df_box", 
                 ylims = c(6,14), dm_ylims = c(-3.7,4), dm_ybreaks=seq(-3,3,3))





#-------------------------------------------------------------------------------

# Relative Scale

#-------------------------------------------------------------------------------

fig_num = "1" 
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/F",fig_num, "/",sep="")


dm_ylims = c(-0.4,0.4)
dm_ybreaks = seq(-0.4, 0.4, 0.2)


# relative mu
#-------------------------------------------------------------------------------
rand_seed <- 1
na1=5; nb1=5; na2=5; nb2=5;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 10
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(2^2/2)/35 + 0.5
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(2^2/2)/35 + 0.6
# Plot experiments and Difference of Means
ggplots <- plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="rmu_box",data_scale = "relative")
gg1 <- ggplots$gg1 + geom_blank(data=tibble(group=c('A','A'),y=c(8,12)))
save_plot(paste(fig_path, "F",fig_num, "_", "rmu_box","_E1", ".tiff", sep = ""),
          gg1, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])

# relative sigma
#-------------------------------------------------------------------------------
# Normalize points
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1^2/2)/0.25e5 + 2e-4
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1^2/2)/0.25e5 + 2e-4
# Plot experiments and Difference of Means
ggplots <- plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="rsigma_box", data_scale = "relative")
gg2 <- ggplots$gg2 + scale_y_continuous(labels=function(x)
  {str_replace(formatC(x, format = "e", digits = 1), 'e([-+])0','e\\1')})
save_plot(paste(fig_path, "F", fig_num, "_", "rsigma_box","_E2", ".tiff", sep = ""),
          gg2, dpi = 600, base_height = ggsize[1], base_width = ggsize[2]*1.1)

# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na3=35; nb3=35;
a3 <- rnorm(na3, mean = 0, sd = 1)
b3 <- rnorm(nb3, mean = 0, sd = 1)
# Normalize points
xa3 <- (a3-mean(a3))/sd(a3)*sqrt(10.5^2/2) + 100
xb3 <- (b3-mean(b3))/sd(b3)*sqrt(10.5^2/2) + 100
# Plot experiments and Difference of Means
plot_experiments(xa1, xb1, xa3, xb3, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="df_box", data_scale = "relative")







