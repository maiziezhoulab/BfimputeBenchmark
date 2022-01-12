#' Evaluation using Monocle
#'
#' @param ematrix Expression single-cell matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param cell_type A vector of cell types with each element representing a
#' cell type for the cell in the corresponding position in \code{ematrix}
#' @param name Title of the plot.
#'
#' @return a list containing the plot of trajectory, the plot of trajectory
#' pseudotime, POS score, and Kendall’s rank correlation score
#' @export
#' @import ggplot2 cowplot ggthemes ggpubr gridExtra ggbeeswarm monocle TSCAN
#'
#' @examples
#'
Monocle_trajectory <- function(ematrix, cell_type, name){
  p.df <- data.frame(CellType = cell_type,
                     row.names = colnames(ematrix))
  pd <- new('AnnotatedDataFrame', data = p.df)
  f.df <- data.frame(gene_short_name = rownames(ematrix),
                     row.names = rownames(ematrix))
  fd <- new('AnnotatedDataFrame', data = f.df)

  tCDS <- newCellDataSet(as.matrix(ematrix),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

  tCDS <- detectGenes(tCDS, min_expr = 0.1)
  # print(head(fData(tCDS)))

  tCDS <- estimateSizeFactors(tCDS)
  # tCDS <- estimateDispersions(tCDS)

  expressed_genes <- rownames(subset(fData(tCDS), num_cells_expressed >= 50))
  # print(head(pData(tCDS)))

  diff_test_res <- differentialGeneTest(tCDS[expressed_genes,],
                                        fullModelFormulaStr = "~CellType")
  ordering_genes <- rownames(subset(diff_test_res, qval < 0.01))

  tCDS <- setOrderingFilter(tCDS, ordering_genes)

  tCDS <- reduceDimension(tCDS, max_components = 2, reduction_method = 'ICA')

  # tCDS <- orderCells(tCDS)
  tCDS <- orderCells(tCDS, num_paths=2, reverse=TRUE)

  tCDS.Pseudotime = pData(tCDS)$Pseudotime

  order.monocle = rank(tCDS.Pseudotime)

  kendall.score = cor(cbind(cell_type.num,order.monocle), method="kendall")[1,2]
  if(kendall.score < 0 ){
    order.monocle = max(order.monocle) + 1 - order.monocle
    kendall.score = cor(cbind(cell_type.num,order.monocle), method="kendall")[1,2]
  }
  subpopulation <- data.frame(cell = colnames(ematrix), sub = cell_type.num)
  orders = list(order(order.monocle))
  POS.score = TSCAN::orderscore(subpopulation, orders)

  trajectoryMonocle = plot_spanning_tree(tCDS, color_by = "CellType",
                                              show_backbone=T,
                                              backbone_color="black") +
    ggtitle(paste0(name,
                   " (POS: ",round(POS.score,3),
                   ", Cor: ", round(kendall.score,3), ")"))

  pseudotimeMonocle = plot_Pseudotime(order.monocle, name, cell_type)
  return(list(trajectoryMonocle = trajectoryMonocle,
              pseudotimeMonocle = pseudotimeMonocle,
              POS.score = POS.score,
              kendall.score = kendall.score))
}

plot_spanning_tree <- function(cds,
                               x=1,
                               y=2,
                               color_by="CellType",
                               show_tree=TRUE,
                               show_backbone=TRUE,
                               backbone_color="black",
                               markers=NULL,
                               show_cell_names=FALSE,
                               cell_name_size=1){
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)

  #print (lib_info_with_pseudo)
  S_matrix <- reducedDimS(cds)

  if (is.null(S_matrix)){
    stop("You must first call reduceDimension() before using this function")
  }

  ica_space_df <- data.frame(t(S_matrix[c(x,y),]))
  colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- minSpanningTree(cds)

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }


  edge_list <- as.data.frame(igraph::get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")

  edge_df <- merge(ica_space_with_state_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)

  edge_df <- plyr::rename(edge_df, c("ICA_dim_1"="source_ICA_dim_1", "ICA_dim_2"="source_ICA_dim_2"))
  edge_df <- merge(edge_df, ica_space_with_state_df[,c("sample_name", "ICA_dim_1", "ICA_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("ICA_dim_1"="target_ICA_dim_1", "ICA_dim_2"="target_ICA_dim_2"))
  edge_df = edge_df[which(!is.na(edge_df$CellType)),]#
  diam <- as.data.frame(as.vector(igraph::V(dp_mst)[igraph::get.diameter(dp_mst, weights=NA)]$name))
  colnames(diam) <- c("sample_name")
  diam <- plyr::arrange(merge(ica_space_with_state_df,diam, by.x="sample_name", by.y="sample_name"), Pseudotime)

  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- melt(exprs(cds[row.names(markers_fData),]))
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "Var1", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    edge_df <- merge(edge_df, markers_exprs, by.x="sample_name", by.y="Var2")
    #print (head(edge_df))
    g <- ggplot(data=edge_df, aes(x=source_ICA_dim_1, y=source_ICA_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=edge_df, aes(x=source_ICA_dim_1, y=source_ICA_dim_2))
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(xend="target_ICA_dim_1", yend="target_ICA_dim_2", color=color_by), size=.3, linetype="solid", na.rm=TRUE)
  }

  g <- g + geom_point(aes_string(color=color_by), size=I(0.5), na.rm=TRUE)
  if (show_backbone){
    #print (diam)
    g <- g + geom_path(aes(x=ICA_dim_1, y=ICA_dim_2), color=I(backbone_color), size=0.75, data=diam, na.rm=TRUE) +
      geom_point(aes_string(x="ICA_dim_1", y="ICA_dim_2", color=color_by), size=I(0.5), data=diam, na.rm=TRUE)
  }

  if (show_cell_names){
    g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
  }
  g <- g +
    # #scale_color_brewer(palette="Set1") +
    # theme(panel.border = element_blank(), axis.line = element_line()) +
    # theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    # theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    # ylab("Component 1") + xlab("Component 2") +
    # theme(legend.position="top", legend.key.height=unit(0.35, "in")) +
    # #guides(color = guide_legend(label.position = "top")) +
    # theme(legend.key = element_blank()) +
    # theme(panel.background = element_rect(fill='white')) +
    scale_color_tableau() +
    # ggtitle(name) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 14),
          legend.title=element_blank())
  return(g)
}

plot_Pseudotime <- function(order.monocle, name, cell_type){
  df = data.frame(order.monocle, cell_type)

  pTSCAN = ggplot(df, aes(x = order.monocle, y = cell_type, colour = cell_type)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() +
    ggtitle(name) +
    theme_classic() +
    xlab("Monocle pseudotime") + ylab("Timepoint") +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 12),
          legend.title=element_blank())

  return(pTSCAN)
}


