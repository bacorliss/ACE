
# Load required packages
#-------------------------------------------------------------------------------
if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(ggplot2)
p_load(cowplot)
library(forcats)

fig_num = "4"
base_dir = "mdm_z"
fig_path = file.path(getwd(), paste(base_dir, "/figure/F",fig_num,sep=""))
dir.create(fig_path, showWarnings = FALSE,recursive = TRUE)


# 7840808, T4
best_df = data.frame(group = c("Progression","BMY22089"), mean = c(0.713, 0.617), 
                     sd = c(0.297/sqrt(8), 0.45/sqrt(8)), n = c(8, 8))
best_df$group <- fct_inorder(best_df$group)
# levels(best_df$group) <- c("Progression Cotrol","BMY22089")
gg<- ggplot(best_df, aes(x = group, y = mean)) +
  geom_bar(stat='identity', fill='white', color="black", width = .5)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=0.2)+
  ylab(expression(paste("Plaque Area (", mm^2,")"))) + xlab("") +
  expand_limits(y = c(0, max(best_df$mean)+1.2*max(best_df$sd))) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=n, y=0.2*max(mean)), size=3) +
  theme_classic(base_size=8)
print(gg)
save_plot(paste(fig_path, "/", fig_num, "1_best.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) 





#11790777, T1
mid_df = data.frame(group = c("FxrWT","Fxr-/-"), mean = c(23.1, 19.4), 
                     sd = c(11.3/sqrt(7), 7.2/sqrt(8)), n = c(7, 8))
mid_df$group <- fct_inorder(mid_df$group)
gg<- ggplot(mid_df, aes(x=group, y=mean)) +
  geom_bar(stat='identity', fill='white', color="black", width = .5)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=0.2)+
  ylab(expression(paste("Plaque Area (%)"))) + xlab("") +
  expand_limits(y = c(0, max(mid_df$mean)+1.2*max(mid_df$sd))) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=n, y=0.2*max(mean)), size=3) +
  # scale_x_discrete(labels = c(
  #   bquote(atop("LDLR"^{"-/-"},"Fxr"^{"WT"})),
  #   bquote(atop("LDLR"^{"-/-"}," Fxr"^{"-/-"})))) +
  theme_classic(base_size=8)
print(gg)
save_plot(paste(fig_path, "/", fig_num, "2_mid.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) 
 
#27683551, F5E
worst_df = data.frame(group = c("-Abx","+Abx"), mean = c(9956, 17196), 
                    sd = c(11578/sqrt(20), 13373/sqrt(18)), n = c(20, 18))
worst_df$group <- fct_inorder(worst_df$group)
gg<- ggplot(worst_df, aes(x=group, y=mean)) +
  geom_bar(stat='identity', fill='white', color="black", width = .5)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=0.2)+
  ylab(expression(paste("Plaque Area (", mu*m^2,")"))) + xlab("") +
  expand_limits(y = c(0, max(worst_df$mean)+1.2*max(worst_df$sd))) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=n, y=0.2*max(mean)), size=3) +
  # scale_x_discrete(labels = c(
  #   bquote(atop("ApoE"^{"-/-"},"Chow~-Abx")),
  #   bquote(atop("ApoE"^{"-/-"},"Chow~+Abx")))) +
  theme_classic(base_size=8)
print(gg)
save_plot(paste(fig_path, "/", fig_num, "3_worst.tiff",sep=""),
          gg, ncol = 1, nrow = 1, base_height = 1.5,
          base_asp = 3, base_width = 1.5, dpi = 600) 
















