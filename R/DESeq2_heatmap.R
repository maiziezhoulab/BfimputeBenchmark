#' Heatmap of some seleted genes of bulk and single-cell expression.
#'
#' Heatmap of some seleted genes of bulk and single-cell expression. The
#' selected genes are Differential Expression Genes detected from bulk using
#' DESeq2
#'
#' @param bulk Expression bulk matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param bulk_cell_type A vector of cell types with each element representing
#' a cell type for the cell in the corresponding position in \code{bulk}
#' @param ematrix Expression single-cell matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param sc_cell_type A vector of cell types with each element representing a
#' cell type for the cell in the corresponding position in \code{ematrix}
#' @param name Title of the plot.
#' @param cell_list A two-element vector containing cell types which are waiting
#' for analysis using DESeq2
#'
#' @return a list containing the heatmap of bulk and single-cell
#' @export
#' @import SingleCellExperiment DESeq2 pheatmap gridExtra cowplot ggplot2
#'
#' @examples
#'
DESeq2_heatmap <- function(bulk, bulk_cell_type,
                           ematrix, sc_cell_type, name,
                           cell_list = c("DEC","H1")) {
  bulk_name = colnames(bulk)
  bulk_keep = sapply(seq(bulk_name), function(i){
    if(bulk_cell_type[i] %in% cell_list) return(TRUE)
    else return(FALSE)
  })
  bulk_de = t(subset(t(bulk),bulk_keep))
  bulk_name = subset(bulk_name,bulk_keep)
  bulk_name_split = subset(bulk_cell_type,bulk_keep)

  bulk_fact = factor(bulk_name_split)
  bulk_deseq = DESeqDataSetFromMatrix(round(bulk_de), DataFrame(bulk_fact), ~ bulk_fact)
  bulk_deseq = DESeq(bulk_deseq)
  bulk_res = results(bulk_deseq)

  #####-------------------------changeable-------------------------#####
  de_remain_subset = bulk_res$padj < 0.05 & (bulk_res$log2FoldChange > 1 | bulk_res$log2FoldChange < -1)
  fold_sub <- subset(bulk_res$log2FoldChange,de_remain_subset)
  de_remain_index_200 = order(abs(fold_sub), decreasing = TRUE)[1:200]

  #####-------------------------draw-------------------------#####
  picture_list = list()
  annocol = data.frame(bulk_fact)
  colnames(annocol) = "Cell type"
  bulk_de_draw <- log10(subset(bulk_de,de_remain_subset)[de_remain_index_200,]+1)
  gap = 0.05
  bks = seq(0, ceiling(max(bulk_de_draw)/gap)*gap, by = gap)
  picture_list$p_bulk = pheatmap(bulk_de_draw, cluster_rows=TRUE, show_rownames=FALSE, treeheight_row = 0,
                                 cluster_cols=FALSE, show_colnames=FALSE, annotation_names_col = FALSE,
                                 annotation_col=as.data.frame(annocol),
                                 breaks = bks,
                                 main = "bulk")

  #####-------------------------sc-------------------------#####
  sc_name = colnames(ematrix)
  sc_keep = sapply(seq(ematrix), function(i){
    if(sc_cell_type[i] %in% cell_list) return(TRUE)
    else return(FALSE)
  })
  sc_name = subset(sc_name,sc_keep)
  sc_name_split = subset(sc_cell_type,sc_keep)

  annocol_sc = as.data.frame(factor(sc_name_split))
  colnames(annocol_sc) = "Cell type"

  ##--##
  sc_de_draw <- log10(subset(t(subset(t(ematrix),sc_keep)),de_remain_subset)[de_remain_index_200,]+1)[picture_list$p_bulk$tree_row$order,]

  picture_list$p_raw = pheatmap(sc_de_draw, cluster_rows=FALSE, show_rownames=FALSE,
                                cluster_cols=FALSE, show_colnames=FALSE, annotation_names_col = FALSE,
                                annotation_col=as.data.frame(annocol_sc),
                                breaks = bks,
                                main = name)
  return(picture_list)
}


#####--------------------------------------------------#####
