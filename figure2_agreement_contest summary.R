



# OVerall contest heatmaps
#------------------------------------------------------------------------

# Re-make heatmap with rectangles based on the selection
my_palette <- colorRampPalette(c(rgb(255, 0, 0,maxColorValue = 255),"white",
                                 rgb(47, 117, 181,maxColorValue = 255)))(n = 299)
col_breaks = c(seq(-1, -.1, length=100), seq(-.09, 0.09, length=100), 
               seq(0.1, 1.0,length=100))
# Function for making selection rectangles around selection cells
makeRects <- function(cells,lwd){
  coords = expand.grid(dim(cells)[1]:1, 1:dim(cells)[2])[cells,]
  xl=coords[,2]-0.49; yb=coords[,1]-0.49; xr=coords[,2]+0.49; yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}
make2RectGroups <- function(cells1,lwd1, cells2, lwd2){
  makeRects(cells1,lwd1)
  makeRects(cells2,lwd2)
}




# Export summary stats for un-scaled data
#-------------------------------------------------------------------------------
dfs_scale <- c(ls(pat="dfs_[1-3]a"),ls(pat="dfs_[1-3]b"))
scale_norm_ind = rep(c(1,3,1,3),2)

# Extract means for each group and subtract from 0.5 (random)
scale_means_from_0.5 <- 0.5 - sapply(dfs_scale, function(x) get(x)$df_plotted$mean)
scale_scores <- t(t(scale_means_from_0.5)/
                    scale_means_from_0.5[cbind(scale_norm_ind,seq_along(scale_norm_ind))])
# Export csv
rownames(scale_scores) <- effect_size_dict[[4]]
# Get statistical significance
scale_scores_sig <- !sapply(dfs_scale, function(x) get(x)$df_plotted$is_mean_0.5) 
scale_score_norm <- sapply(scale_norm_ind, function(ind,len) 
  ifelse(1:len == ind, TRUE,FALSE), length(effect_size_dict[[4]]))

png(paste("figure/F", fig_num, "/F", fig_num, "es_contest scale.png",sep=""),    
    width = 5*300, height = 2.55*300, res = 300, pointsize = 8)  
heatmap.2(scale_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {make2RectGroups(scale_scores_sig,1,scale_score_norm,3)}, 
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(effect_size_dict[[4]], function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(scale_scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(scale_scores)[1]),
          notecol="black",notecex=1, lwid=c(0.001,5),lhei=c(0.001,5),margins =c(0,0))
dev.off()

# Export summary stats for relative scale data
#-------------------------------------------------------------------------------
# Export summary stats for relative scale data
dfs_rscale <- c(ls(pat="dfs_[4-6]a"),ls(pat="dfs_[4-6]b"))
rscale_norm_ind = rep(c(2,4,2,4),2)
# Extract means for each group and subtract from 0.5 (random)
rscale_means_from_0.5 <- 0.5 - sapply(dfs_rscale, function(x) get(x)$df_plotted$mean)
rscale_scores <- t(t(rscale_means_from_0.5)/
                     rscale_means_from_0.5[cbind(rscale_norm_ind,seq_along(rscale_norm_ind))])
# Export csv
rownames(rscale_scores) <- effect_size_dict[[2]]
# Get statistical significance
rscale_scores_sig <- !sapply(dfs_rscale, function(x) get(x)$df_plotted$is_mean_0.5) 
rscale_score_norm <- sapply(rscale_norm_ind, function(ind,len) 
  ifelse(1:len == ind, TRUE,FALSE), length(effect_size_dict[[4]]))

png(paste("figure/F", fig_num, "/F", fig_num, "es_contest relative scale.png",sep=""),    
    width = 5*300, height = 2.55*300, res = 300, pointsize = 8)  
heatmap.2(rscale_scores, trace = "none", dendrogram = "none", key = FALSE,
          add.expr = {make2RectGroups(rscale_scores_sig,1,rscale_score_norm,3)}, 
          col = my_palette,  Rowv=F, Colv=F, sepwidth=c(0,0),
          labRow =  sapply(effect_size_dict[[4]], function(x) parse(text=x)),labCol = "",
          cellnote=matrix(sapply(rscale_scores,function(x) sprintf("%0.2+f",x)),
                          nrow = dim(rscale_scores)[1]),
          notecol="black",notecex=1, lwid=c(0.01,5),lhei=c(0.01,5),margins =c(0,0))
dev.off()

