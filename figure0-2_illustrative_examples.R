 


library(ggplot2)
library(dplyr)
library(cowplot)



# User defined functions
plot_experiments <- function(xa1,xb1,xa2,xb2, fig_num, fig_path, base_name, scale = "absolute") {
  
  blank_df <- tibble(group=c("A","A"), y = c(max(c(xa1,xb1,xa2,xb2)), 
                                             min(c(xa1,xb1,xa2,xb2))))
  # Experiment 1
  df <- tibble(group = as.factor(c(rep("A",na1), rep("B",nb1))), y = c(xa1,xb1))
  df_stats <- df %>% group_by(group) %>% summarise(mean = mean(y),sd = sd(y))
  gg1 <- ggplot(data=df, aes(x=group, y=y)) +
    geom_jitter(width = 0.2, size=0.3, color = "gray50") + 
    geom_point(data = df_stats, aes(x=group, y=mean), pch=3) +
    geom_blank(data=blank_df, aes(x=group, y=y)) +
    theme_classic(base_size=8) + xlab("Group") + ylab("f(x)") +
    geom_errorbar(data = df_stats, aes(x=group,ymin=mean-sd,
                                       ymax=mean+sd, y=mean), width=0.4) +
    expand_limits(y = min(df$y) + 1.2*(max(df$y)-min(df$y)))
  gg1
  save_plot(paste(fig_path, "F",fig_num, "_", base_name,"_E1", ".tiff", sep = ""),
            gg1, dpi = 600, base_height = 1.2, base_width = .8)
  

  # Experiment 2
  df <- tibble(group = as.factor(c(rep("A",na2),rep("B",nb2))),y = c(xa2,xb2))
  df_stats <- df %>% group_by(group) %>% summarise(mean = mean(y),sd = sd(y))
  gg2 <- ggplot(data=df, aes(x=group, y=y)) +
    geom_jitter(width = 0.2, size=0.3, color = "gray50") + 
    geom_point(data = df_stats, aes(x=group, y=mean), pch=3) +
    geom_blank(data=blank_df, aes(x=group, y=y)) +
    theme_classic(base_size=8) + xlab("Group") + ylab("f(x)") +
    geom_errorbar(data = df_stats, aes(x=group,ymin=mean-sd,
                                       ymax=mean+sd, y=mean), width=0.4) +
    expand_limits(y = min(df$y) + 1.2*(max(df$y)-min(df$y)))
  gg2
  save_plot(paste(fig_path, "F", fig_num, "_", base_name,"_E2", ".tiff", sep = ""),
            gg2, dpi = 600, base_height = 1.2, base_width = .8)
  
  
  # Difference in Means
  
  # Stats for each group
  dfo <- tibble(mu_d1 = mean(xb1) - mean(xa1),sigma_d1 = sqrt(sd(xa1)^2 + sd(xb1)^2),
                mu_dm1 = mean(xb1) - mean(xa1), sigma_dm1 = sqrt(sd(xa1)^2/na1 + sd(xb1)^2/nb1),
                mu_d2 = mean(xb2) - mean(xa2), sigma_d2 = sqrt(sd(xa2)^2 + sd(xb2)^2),
                mu_dm2 = mean(xb2) - mean(xa2), sigma_dm2 = sqrt(sd(xa2)^2/na2 + sd(xb2)^2/nb2),
                rmu_dm1 = mu_dm1/mean(xa1), rsigma_dm1 = sigma_dm1/mean(xa1),
                rmu_dm2 =    mu_dm2/mean(xa2), rsigma_dm2 = sigma_dm2/mean(xa2),
                df1=length(xa1)+length(xb1)-2, df2=length(xa2)+length(xb2)-2)

  # n in this case is number of samples for distribution histogram
  n1 = 1e4; n2 = 1e4
  if (scale=="relative") {
    df_sample <- tibble(group = as.factor(c(rep(1,n1),rep(2,n2))), 
                        y = c(rnorm(n1,mean = dfo$mu_dm1, sd = dfo$sigma_dm1)/mean(xa1),
                              rnorm(n2,mean = dfo$mu_dm2, sd = dfo$sigma_dm2)/mean(xa2)))
    # Both experiments
    gg3 <- ggplot(data = df_sample, aes(x=group, y=y)) +
      geom_violin(fill = "gray95") + theme_classic(base_size=8) + 
      xlab("Exp") + ylab(expression(mu[DM]*phantom(.) / phantom(.) * mu[A]))
    # Exp 1 only
    gg4 <- ggplot(data = subset(df_sample, group==1), aes(x=group, y=y)) +
      geom_violin(fill = "gray95") + theme_classic(base_size=8) + 
      xlab("Exp") + ylab(expression(mu[DM]))
    gg5 <- ggplot(data = subset(df_sample, group==2), aes(x=group, y=y)) +
      geom_violin(fill = "gray95") + theme_classic(base_size=8) + 
      xlab("DM") + ylab(NULL) +
      theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
            plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
  } else {y_label = expression(mu[DM]); denominator = 1; 
  df_sample <- tibble(group = as.factor(c(rep(1,n1),rep(2,n2))), 
                      y = c(rnorm(n1,mean = dfo$mu_dm1, sd = dfo$sigma_dm1)/1,
                            rnorm(n2,mean = dfo$mu_dm2, sd = dfo$sigma_dm2)/1))
  # Both experiments
  gg3 <- ggplot(data = df_sample, aes(x=group, y=y)) +
    geom_violin(fill = "gray95") + theme_classic(base_size=8) + 
    xlab("Exp") + ylab(expression(mu[DM]))
 # Exp 1 only
  gg4 <- ggplot(data = subset(df_sample, group==1), aes(x=group, y=y)) +
    geom_violin(fill = "gray95") + theme_classic(base_size=8) + 
    xlab("DM") + ylab(NULL) +
    theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
  gg5 <- ggplot(data = subset(df_sample, group==2), aes(x=group, y=y)) +
    geom_violin(fill = "gray95") + theme_classic(base_size=8) + 
    xlab("DM") + ylab(NULL) +
    theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
  }
  gg3
  save_plot(paste(fig_path, "F", fig_num, "_", base_name,"_D", ".tiff", sep = ""),
            gg3, dpi = 600, base_height = 1.1, base_width = 1)
  save_plot(paste(fig_path, "F", fig_num, "_", base_name,"_Exp1", ".tiff", sep = ""),
            gg4, dpi = 600, base_height = 1, base_width = .4)
  save_plot(paste(fig_path, "F", fig_num, "_", base_name,"_Exp2", ".tiff", sep = ""),
            gg5, dpi = 600, base_height = 1, base_width = .4)
  print(t(dfo))
}


#-------------------------------------------------------------------------------

# Single Result

#-------------------------------------------------------------------------------
fig_num = "1" 
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/F",fig_num, "/",sep="")


# Illustrative group
rand_seed <- 2
na1=6; nb1=6; na2=6; nb2=6;
rand_seed <- 0
na1=6; nb1=6; na2=6; nb2=6;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) * 1.4 + 0
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) * 1.4 + 1
xa2 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) * 1 + 0
xb2 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) * 0.5 + 0
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2,fig_num=fig_num, fig_path=fig_path, base_name="illustrative_box")

# 
# # Mean offset
# rand_seed <- 0
# na1=6; nb1=6; na2=6; nb2=6;
# rand_seed <- 0
# na1=6; nb1=6; na2=6; nb2=6;
# a1 <- rnorm(na1, mean = 0, sd = 1)
# b1 <- rnorm(nb1, mean = 0, sd = 1)
# a2 <- rnorm(na2, mean = 0, sd = 1)
# b2 <- rnorm(nb2, mean = 0, sd = 1)
# # Normalize points
# xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2)
# xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 1
# xa2 <- (a2-mean(a2))/sd(a2)
# xb2 <- (b2-mean(b2))/sd(b2)
# # Plot experiments and Difference of Means
# plot_experiments(xa1,xb1,xa2,xb2,fig_num=fig_num, fig_path=fig_path, base_name="mu_single_box")
# 
# 
# 
# # std
# rand_seed <- 0
# na1=6; nb1=6; na2=6; nb2=6;
# rand_seed <- 0
# na1=6; nb1=6; na2=6; nb2=6;
# a1 <- rnorm(na1, mean = 0, sd = 1)
# b1 <- rnorm(nb1, mean = 0, sd = 1)
# a2 <- rnorm(na2, mean = 0, sd = 1)
# b2 <- rnorm(nb2, mean = 0, sd = 1)
# # Normalize points
# xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) + 0
# xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 0
# xa2 <- (a2-mean(a2))/sd(a2)
# xb2 <- (b2-mean(b2))/sd(b2)
# # Plot experiments and Difference of Means
# plot_experiments(xa1,xb1,xa2,xb2,fig_num=fig_num, fig_path=fig_path, base_name="std_single_box")




#-------------------------------------------------------------------------------

# Absolute Scale

#-------------------------------------------------------------------------------

fig_num = "8" 
dir.create(file.path(getwd(), paste("figure/SF",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/SF",fig_num, "/",sep="")

# Decreased mean offset for experiment 2
rand_seed <- 0
na1=6; nb1=6; na2=6; nb2=6;
# Decreased mean offset for experiment 2
rand_seed <- 0
na1=6; nb1=6; na2=6; nb2=6;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 14
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1^2/2) + 12
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2,fig_num=fig_num, fig_path=fig_path, base_name="mu_box")





# Decreased standard deviation for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na1=6; nb1=6; na2=6; nb2=6;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(2^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(2^2/2) + 10
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1^2/2) + 10
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, base_name="sigma_box")



# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na1=6; nb1=6; na2=21; nb2=21;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 10
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1^2/2) + 10
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, base_name="df_box")






#-------------------------------------------------------------------------------

# Relative Scale

#-------------------------------------------------------------------------------

fig_num = "1" 
dir.create(file.path(getwd(), paste("figure/F",fig_num,sep="")), showWarnings = FALSE)
fig_path = paste("figure/F",fig_num, "/",sep="")


# relative mu
#-------------------------------------------------------------------------------
rand_seed <- 0
na1=6; nb1=6; na2=6; nb2=6;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 15
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(2^2/2) + 20
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(2^2/2) + 25
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, base_name="rmu_box",scale = "relative")


# relative sigma
#-------------------------------------------------------------------------------
rand_seed <- 0
na1=6; nb1=6; na2=6; nb2=6;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(2^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(2^2/2) + 10
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(2^2/2) + 20
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(2^2/2) + 20
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, base_name="rsigma_box", scale = "relative")



# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na1=6; nb1=6; na2=21; nb2=21;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1^2/2) + 10
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(2^2/2) + 20
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(2^2/2) + 20
# Plot experiments and Difference of Means
plot_experiments(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, base_name="df_box",scale = "relative")







