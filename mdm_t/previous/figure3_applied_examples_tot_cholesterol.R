
# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(ggplot2)
p_load(cowplot)

fig_num = "3"
base_dir = "mdm_t"
fig_path = file.path(getwd(), paste(base_dir, "/figure/F",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE,recursive = TRUE)


# 24797634, 1C
best_df = data.frame(group = c("Alfp-cre","Alfp-creTRBfl/fl"), mean = c(3.45, 3.26), 
                     sd = c(0.24/sqrt(6), 0.22/sqrt(6)), n = c(6, 6))
gg<- ggplot(best_df, aes(x=group, y=mean)) +
  geom_bar(stat='identity', fill='white', color="black", width = .5)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=0.2)+
  ylab("Total Cholesterol (mmol/L)   ") + xlab("") +
  expand_limits(y = c(0, max(best_df$mean)+1.2*max(best_df$sd))) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=n, y=0.2*max(mean)), size=3) +
  # scale_x_discrete(labels = c(bquote("Alfp-cre"),bquote("Alfp-cre"~"TR"*beta^{fl/fl}))) +
  theme_classic(base_size=8)
print(gg)
save_plot(paste(fig_path, "/", fig_num, "1_best.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) 





#11790777, T1
mid_df = data.frame(group = c("ApoE-/-CD6WT","ApoE-/-CD6-/-"), mean = c(2518, 2876), 
                     sd = c(257/sqrt(8), 506/sqrt(6)), n = c(8, 6))
gg<- ggplot(mid_df, aes(x=group, y=mean)) +
  geom_bar(stat='identity', fill='white', color="black", width = .5)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=0.2)+
  ylab("Total Cholesterol (mg/dL)  ") + xlab("") +
  expand_limits(y = c(0, max(mid_df$mean)+1.2*max(mid_df$sd))) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=n, y=0.2*max(mean)), size=3) +
  # scale_x_discrete(labels = c(bquote("ApoE"^{"-/-"}),bquote("ApoE"^{"-/-"}*"CD6"^{"-/-"}))) +
  theme_classic(base_size=8)
print(gg)
save_plot(paste(fig_path, "/", fig_num, "2_mid.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) 
 
#27683551, F5E
worst_df = data.frame(group = c("ApoE-/-","ApoE-/- + PAO"), mean = c(568, 742), 
                    sd = c(81/sqrt(6), 256/sqrt(4)), n = c(6, 4))
gg<- ggplot(worst_df, aes(x=group, y=mean)) +
  geom_bar(stat='identity', fill='white', color="black", width = .5)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=0.2)+
  ylab("Total Cholesterol (mg/dL)  ") + xlab("") +
  expand_limits(y = c(0, max(worst_df$mean)+1.2*max(worst_df$sd))) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=n, y=0.2*max(mean)), size=3) +
  # scale_x_discrete(labels = c(bquote("ApoE"^{"-/-"}),bquote("ApoE"^{"-/-"}~"+PAO"))) +
  theme_classic(base_size=8)
print(gg)
save_plot(paste(fig_path, "/", fig_num, "3_worst.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) 
















