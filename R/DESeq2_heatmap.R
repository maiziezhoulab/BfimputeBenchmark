library(SingleCellExperiment)
library(DESeq2)
#library(MAST)
library(pheatmap)
library(gridExtra)
library(cowplot)
library(ggplot2)

bulk_name = colnames(bulk)
bulk_keep = sapply(1:dim(bulk)[2], function(i){
  split_name <- unlist(strsplit(bulk_name[i],"_"))
  if(split_name[1] %in% cell_list) return(TRUE)
  else return(FALSE)
})
bulk_de = t(subset(t(bulk),bulk_keep))
bulk_name = subset(bulk_name,bulk_keep)
bulk_name_split = sapply(bulk_name, function(i){
  split_name <- unlist(strsplit(i,"_"))[1]
})

bulk_fact = factor(bulk_name_split)
bulk_deseq = DESeqDataSetFromMatrix(round(bulk_de), DataFrame(bulk_fact), ~ bulk_fact)
bulk_deseq = DESeq(bulk_deseq)
bulk_res = results(bulk_deseq)

#####-------------------------changeable-------------------------#####
de_remain_subset = bulk_res$padj < 0.05 & (bulk_res$log2FoldChange > 1 | bulk_res$log2FoldChange < -1)
fold_sub <- subset(bulk_res$log2FoldChange,de_remain_subset)
de_remain_index_200 = order(abs(fold_sub), decreasing = TRUE)[1:200]

#####-------------------------draw-------------------------#####
annocol = data.frame(bulk_fact)
colnames(annocol) = "Cell type"
bulk_de_draw <- log10(subset(bulk_de,de_remain_subset)[de_remain_index_200,]+1)
bks = seq(0, ceiling(max(bulk_de_draw)/gap)*gap, by = gap)
picture_list$p_bulk = pheatmap(bulk_de_draw, cluster_rows=TRUE, show_rownames=FALSE, treeheight_row = 0,
                               cluster_cols=FALSE, show_colnames=FALSE, annotation_names_col = FALSE,
                               annotation_col=as.data.frame(annocol),
                               breaks = bks,
                               main = "bulk")

#####-------------------------sc-------------------------#####
sc_name = colnames(res)
sc_keep = sapply(1:ncol(res), function(i){
  split_name <- unlist(strsplit(sc_name[i],"_"))
  if(split_name[1] %in% cell_list) return(TRUE)
  else return(FALSE)
})
sc_name = subset(sc_name,sc_keep)
sc_name_split = sapply(sc_name, function(i){
  split_name <- unlist(strsplit(i,"_"))[1]
})

annocol_sc = as.data.frame(factor(sc_name_split))
colnames(annocol_sc) = "Cell type"

##--##
sc_de_draw <- log10(subset(t(subset(t(res),sc_keep)),de_remain_subset)[de_remain_index_200,]+1)[picture_list$p_bulk$tree_row$order,]

picture_list$p_raw = pheatmap(sc_de_draw, cluster_rows=FALSE, show_rownames=FALSE,
                              cluster_cols=FALSE, show_colnames=FALSE, annotation_names_col = FALSE,
                              annotation_col=as.data.frame(annocol_sc),
                              breaks = bks,
                              main = name)


#####--------------------------------------------------#####
