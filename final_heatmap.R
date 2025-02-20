library (pheatmap)
library (gplots)
library(gridExtra)

library(grid)

median_mat <- readRDS("/home/xkdl27/mouse_DB/heat_map/median_expression_matrix.rds")

TRIM <- read.table('mouse_DB/cancer_expression_DB/Mouse_trim_gene_list.txt', sep = '\t', header = T)

breaksList = seq(-1.5, 1.5, by = 0.1)
breaksList=append(breaksList, 2)
breaksList=append(breaksList, -2, 0)
mycol <- colorpanel(n=length(breaksList)-1,low="blue", mid = "white", high="red")


normalized_Trim_matrix <- median_mat[rownames(median_mat) %in% TRIM$Symbol,] 

normal_matrix <- normalized_Trim_matrix[, 1:44]
cancer_matrix <- normalized_Trim_matrix[, 45:65]


res1 = pheatmap(log2(normal_matrix+1), show_rownames=F, show_colnames=T, 
                cluster_cols=T, cluster_rows=T, scale="row", col=mycol, breaks = breaksList, 
                clustering_distance_rows ="euclidean", cex=0.8, clustering_method="ward.D", 
                border_color=NA, legend=FALSE)

cancer_matrix_reorder <- cancer_matrix[res1$tree_row$order,]
res2 = pheatmap(log2(cancer_matrix_reorder + 1), show_rownames=T, show_colnames=T, cluster_cols=T, cluster_rows=F, scale="row",
                col=mycol, breaks=breaksList, clustering_distance_rows="euclidean", cex=0.8, clustering_method="ward.D", border_color=NA)


# align res1 & res2

max_height <- unit.pmax(res1$gtable$heights, res2$gtable$heights)
res1$gtable$heights <- max_height
res2$gtable$heights <- max_height


# extract only legend
legend_grob <- res1$gtable$grobs[[5]]



grid.arrange(legend_grob, res1$gtable, res2$gtable, ncol = 3, widths=c(0.15, 1, 1))


grid.newpage()
pushViewport((viewport(x=0.015, y=0.45, just="left")))
grid.draw(legend_grob)


empty_grob <- rectGrob(gp=gpar(col="white"))

grid.arrange(empty_grob, res1$gtable, res2$gtable, ncol = 3, widths=c(0.1, 1, 1))