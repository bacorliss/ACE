 

#' A series of examples to illustrate the three parameters of unscaled agreement
#' and relative agreement. Each example compares the difference in means or 
#' relative difference in means between a hypothetical control and experiment
#' group.

# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
source("R/illustrative_plot.R")

# Figure parameters
#-------------------------------------------------------------------------------
base_dir = "mdm_t"
fig_num = "2_1" 
dir.create(file.path(getwd(), base_dir,"figure"), showWarnings = FALSE)
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(file.path(getwd(), fig_path), showWarnings = FALSE)
grp_color = brewer.pal(3, "Set1")
yaxis_font_size = 7
base_font_size = 8

# Figure sizes, height x length
ggsize = c(1.2, .65)



#-------------------------------------------------------------------------------

# Absolute Scale

#-------------------------------------------------------------------------------

# Decreased mean offset for experiment 2
set.seed(0)
na1=15; nb1=15; na2=15; nb2=15;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)

# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(2.1^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(2.1^2/2) + 10.3
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(2.1^2/2) + 10
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(2.7^2/2) + 11.2
# Plot experiments and Difference of Means
gg_mu <- plot_experiment_pair(xa1, xb1, xa2, xb2, fig_num=fig_num, fig_path=fig_path, base_name="_1_unsc_mu", 
                 ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2))


# Increase alpha
#------------------------------------------------------------------------------
# Keep group 1, change alpha only
gg_mu_alpha <- plot_experiment_pair(xa1, xb1, xa1, xb1, fig_num=fig_num, fig_path=fig_path, base_name="_4_unsc_alpha", 
                                ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2),
                                conf.level.2 = 1-0.05/100)
# Steal jitter positions from mu_plot1 to use for alpha_plot2 to show they are 
# exactly the same points
q1 <- ggplot_build(gg_mu$gg1)
q1$data[[1]]$x
q2 <- ggplot_build(gg_mu_alpha$gg2)
q2$data[[1]]$x <- q1$data[[1]]$x
gq <- ggplot_gtable(q2)
save_plot(paste(fig_path, "/F",fig_num, "_", "_4_unsc_alpha","_E2", ".tiff", sep = ""),
          gq, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])


# Increased std experiment 2
#-------------------------------------------------------------------------------
xa3 <- (a2-mean(a2))/sd(a2)*sqrt(4.^2/2) + 10
xb3 <- (b2-mean(b2))/sd(b2)*sqrt(4.^2/2) + 10.3
# Plot experiments and Difference of Means
plot_experiment_pair(xa1,xb1,xa3,xb3, fig_num=fig_num, fig_path=fig_path, base_name="_2_unsc_sigma", 
                 ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2))


# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na2=6; nb2=6;
a3 <- rnorm(na2, mean = 0, sd = 1)
b3 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa3 <- (a3-mean(a3))/sd(a3)*sqrt(2.1^2/2) + 10
xb3 <- (b3-mean(b3))/sd(b3)*sqrt(2.1^2/2) + 10.3
# Plot experiments and Difference of Means
plot_experiment_pair(xa1,xb1,xa3,xb3, fig_num=fig_num, fig_path=fig_path, base_name="_3_unsc_df", 
                 ylims = c(6,14), dm_ylims = c(-2.25,2.25), dm_ybreaks=seq(-2,2,2))












#-------------------------------------------------------------------------------

# Relative Scale

#-------------------------------------------------------------------------------

dm_ylims = c(-0.2,0.2)
dm_ybreaks = seq(-0.2, 0.2, 0.1)


# Relative mu
#-------------------------------------------------------------------------------
set.seed(0)
na1=15; nb1=15; na2=15; nb2=15;
a1 <- rnorm(na1, mean = 0, sd = 1)
b1 <- rnorm(nb1, mean = 0, sd = 1)
a2 <- rnorm(na2, mean = 0, sd = 1)
b2 <- rnorm(nb2, mean = 0, sd = 1)
# Normalize points
xa1 <- (a1-mean(a1))/sd(a1)*sqrt(1.65^2/2) + 10
xb1 <- (b1-mean(b1))/sd(b1)*sqrt(1.65^2/2) + 10.15
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1.30^2/2)/15 + 0.5
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1.30^2/2)/15 + 0.54
# Plot experiments and Difference of Means
ggs_rmu <- plot_experiment_pair(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="_5_rel_rmu",data_scale = "relative")
gg1 <- ggs_rmu$gg1 + geom_blank(data=tibble(group=c('X','X'),y=c(8,12))) + scale_y_continuous(breaks = seq(8,12,1))
save_plot(paste(fig_path, "/F",fig_num, "_", "_5_rel_rmu","_E1", ".tiff", sep = ""),
          gg1, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])

# Alpha dm resuse group 1 from mu example
# ---------------------
gg_ralpha <- plot_experiment_pair(xa1,xb1,xa1/3.5,xb1/3.5, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="_8_rel_alpha",data_scale = "relative",
                 conf.level.2 = 1-0.05/75)
#Overwrite the figure for group 2 so it has exactly same placement as group 1
q1 <- ggplot_build(gg1)
q2 <- ggplot_build(gg_ralpha$gg2 + coord_cartesian(ylim = c(2.28,3.6)))
q2$data[[1]]$x <- q1$data[[1]]$x
q2$data[[1]]$x[seq(na1+1, na1+na2, 1)] <- q2$data[[1]]$x[seq(na1+1, na1+na2, 1)] - 1
gq <- ggplot_gtable(q2)

# grid.draw(gq)

save_plot(paste(fig_path, "/F",fig_num, "_", "_8_rel_alpha","_E2", ".tiff", sep = ""),
          gq, dpi = 600, base_height = ggsize[1], base_width = ggsize[2])


# Relative sigma
#-------------------------------------------------------------------------------
# Normalize points
xa2 <- (a2-mean(a2))/sd(a2)*sqrt(1^2/2)/0.14e5 + 2.2e-4
xb2 <- (b2-mean(b2))/sd(b2)*sqrt(1^2/2)/0.14e5 + 2.233e-4
# Plot experiments and Difference of Means
ggplots <- plot_experiment_pair(xa1,xb1,xa2,xb2, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="_6_rel_rsigma", data_scale = "relative")
gg2 <- ggplots$gg2 + 
  geom_blank(data=data.frame(x=c(1, 1),y=c(1.8e-4, 3e-4)), aes(x=x, y=y)) +
  scale_y_continuous(breaks = c(1e-4, 2e-4, 3e-4), labels=function(x)
  {str_replace(formatC(x, format = "e", digits = 0), 'e([-+])0','e\\1')})
print(gg2)
save_plot(paste(fig_path, "/F", fig_num, "_", "_6_rel_rsigma","_E2", ".tiff", sep = ""),
          gg2, dpi = 600, base_height = ggsize[1], base_width = ggsize[2]*1.1)


# Increased sample size for experiment 2
#-------------------------------------------------------------------------------
rand_seed <- 0
na3=6; nb3=6;
a3 <- rnorm(na3, mean = 0, sd = 1)
b3 <- rnorm(nb3, mean = 0, sd = 1)
# Normalize points
xa3 <- (a3-mean(a3))/sd(a3)*sqrt(18^2/2) + 100
xb3 <- (b3-mean(b3))/sd(b3)*sqrt(18^2/2) + 101.5
# Plot experiments and Difference of Means
ggplots <- plot_experiment_pair(xa1, xb1, xa3, xb3, fig_num=fig_num, fig_path=fig_path, dm_ylims = dm_ylims,
                 dm_ybreaks = dm_ybreaks, base_name="_7_rel_df", data_scale = "relative")

