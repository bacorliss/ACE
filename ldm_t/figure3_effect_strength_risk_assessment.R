


#' Runs a collection of risk assessments for test how often each candidate statistic 
#' incorrectly predicts which of two results have greater noteworthiness (higher
#' evidence of practical significance)

if (!require("pacman")) {install.packages("pacman")}; library(pacman)
p_load(ggplot2)
p_load(gplots)

# Figure parameters
#-------------------------------------------------------------------------------
fig_num = "2"
base_dir = "ldm_t"
sum_fig_path = paste(getwd(),'/', base_dir,"/figure/F",fig_num, sep="")
dir.create(sum_fig_path, showWarnings = FALSE, recursive = TRUE)

data_path <- paste(getwd(),'/', base_dir,"/temp", sep="")
dir.create(data_path, showWarnings = FALSE, recursive = TRUE)


fig_height_inch <- 2.5
fig_width_inch <- 6.25

# Heatmap formatting
#------------------------------------------------------------------------

# Remake heatmap with rectangles based on the selection
my_palette <- colorRampPalette(c(rgb(47, 117, 181,maxColorValue = 255),"white",
                                 rgb(255, 0, 0,maxColorValue = 255)))(n = 299)
col_breaks = c(seq(-1, -.1, length=100), seq(-.09, 0.09, length=100), 
               seq(0.1, 1.0,length=100))
# Function for making selection rectangles around selection cells
makeRects <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[2])[cells,]
  xl=coords[,2]-0.485; yb=coords[,1]-0.485; xr=coords[,2]+0.48; yt=coords[,1]+0.48
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}

add_underline <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[2])[cells,]
  xl=coords[,2]-0.49; yb=coords[,1]-0.49; xr=coords[,2]+0.49; yt=coords[,1]+0.49
  segments(xl+.07, yb+.16, xr-.08, yb+.16, col = "black", lty = "solid", lwd = lwd)
}



# Input arguments
n_samples = 1E3
use_pseudo_samples = TRUE

# Export summary stats for unscaled data
#-------------------------------------------------------------------------------

if (!file.exists(file.path(data_path, "df_unscaled_pos.RDS"))) {
  print("SFig 9: unscaled agreement contest pos")
  source(file.path(base_dir, "sfigure9_pos_unscaled_risk_assessment_crit.R"))
} else {load(file = file.path(data_path, "df_unscaled_pos.RDS"))}

if (!file.exists(file.path(data_path, "df_unscaled_neg.RDS"))) {
  print("SFig 10: unscaled agreement contest neg")
  source(file.path(base_dir, "sfigure10_neg_unscaled_risk_assessment_crit.R"))
} else {load(file = file.path(data_path, "df_unscaled_neg.RDS"))}


# Heatmap of unscaled data
#-------------------------------------------------------------------------------
dfs_unscaled <- c(df_unscaled_pos,df_unscaled_neg)

# Extract means for each group and subtract from 0.5 (random)
scale_means_from_0.5 <- sapply(dfs_unscaled, function(x) x$df_plotted$mean) - 0.5
scale_scores <- sweep(scale_means_from_0.5, 2, abs(apply(scale_means_from_0.5,2,min)), FUN = '/')

# Export csv
rownames(scale_scores) <- attr(dfs_unscaled[[1]]$df_es, "varnames")
# Identity which cells are statistically significant from random
scale_scores_sig <- !sapply(dfs_unscaled, function(x) x$df_plotted$is_mean_0.5) 

# Zero color to white for fields that are not statistically significant
zeroed_scale_scores <- scale_scores
zeroed_scale_scores[!scale_scores_sig] <- 0
png(paste(sum_fig_path, "/F", fig_num, "_risk_assessment_unscaled_agreement.png",sep=""),    
    width = fig_width_inch*300, height = fig_height_inch*300, res = 300, pointsize = 8)  
heatmap.2(zeroed_scale_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {add_underline(scale_scores_sig,1.5);}, 
          col = my_palette,  Rowv = F, Colv = F, sepwidth = c(0,0),
          labRow =  sapply(attr(dfs_unscaled[[1]]$df_es, "varnames"),
                           function(x) parse(text=x)),labCol = "",
          cellnote = matrix(sapply(scale_scores,function(x) sprintf("%0.2+f",x)),
                            nrow = dim(scale_scores)[1]),
          breaks = col_breaks,
          notecol ="black", notecex = 1, lwid = c(0.001,5),lhei = c(0.001,5),margins = c(0,0))
dev.off()





# Export summary stats for relative scale data
#-------------------------------------------------------------------------------

if (!file.exists(file.path(data_path, "df_relative_pos.RDS"))) {
  print("SFig 10: relative agreement contest pos")
  source(file.path(base_dir, "sfigure11_pos_rel_risk_assessment_crit.R"))
} else {load(file = file.path(data_path, "df_relative_pos.RDS"))}

if (!file.exists(file.path(data_path, "df_relative_neg.RDS"))) {
  print("SFig 11: relative agreement contest neg")
  source(file.path(base_dir, "sfigure12_neg_rel_risk_assessment_crit.R"))
} else {load(file = file.path(data_path, "df_relative_neg.RDS"))}

fig_num = "2" 


# Heatmap of relative data
#-------------------------------------------------------------------------------
# Export summary stats for relative scale data
dfs_relative <- c(df_relative_pos,df_relative_neg)

# Extract means for each group and subtract from 0.5 (random)
relative_means_from_0.5 <- sapply(dfs_relative, function(x) x$df_plotted$mean) - 0.5
# scale_means_from_0.5 <- sapply(dfs_unscaled, function(x) x$df_plotted$mean) - 0.5
rscale_scores <- sweep(relative_means_from_0.5, 2, abs(apply(relative_means_from_0.5,2,min)), FUN = '/')

# Export csv
rownames(rscale_scores) <- attr(dfs_relative[[1]]$df_es,"varnames")
# Get statistical significance
rscale_scores_sig <- !sapply(dfs_relative, function(x) x$df_plotted$is_mean_0.5) 

# Zero color to white for fields that are not statistically significant
zeroed_rscale_scores <- rscale_scores
zeroed_rscale_scores[!rscale_scores_sig] <- 0
png(paste(sum_fig_path,"/F", fig_num, "_risk_assessment_relative_agreement.png",sep=""),    
    width = fig_width_inch*300, height = fig_height_inch*300, res = 300, pointsize = 8)  
heatmap.2(zeroed_rscale_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {add_underline(rscale_scores_sig,1.5);},
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(attr(dfs_relative[[1]]$df_es,"varnames_pretty"),
                           function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(rscale_scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(rscale_scores)[1]),
          breaks = col_breaks,
          notecol="black",notecex=1, lwid=c(0.01,5),lhei=c(0.01,5),margins =c(0,0))
dev.off()

