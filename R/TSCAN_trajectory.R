#' Evaluation using TSCAN
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
#' @import ggplot2 cowplot ggthemes ggpubr gridExtra ggbeeswarm SingleCellExperiment scater TSCAN
#'
#' @examples
#'
TSCAN_trajectory <- function(ematrix, cell_type, name){
  cell.names = colnames(ematrix)
  res <- ematrix
  colnames(res) <- NULL
  rownames(res) <- NULL
  res <- SingleCellExperiment(
    assays = list(logcounts = log2(as.matrix(res) + 1)),
    colData = data.frame(cell_type = cell_type)
  )
  set.seed(123456)
  res <- runPCA(res)

  data.pca = t(reducedDim(res))[1:20,]

  colnames(data.pca) = cell.names
  Kcluster = sum(!is.na(unique(cell_type)))
  K.arrage = Kcluster+1
  data.clust <- TSCAN::exprmclust(data.pca,
                                  reduce = F,
                                  clustermethod = "kmeans",
                                  clusternum = K.arrage)

  data.order <- TSCAN::TSCANorder(data.clust, orderonly = FALSE)

  order.tscan = rep(NA, ncol(ematrix))
  names(order.tscan) = cell.names

  order.tscan[data.order$sample_name] <- data.order$Pseudotime

  cell_type.num.temp = cell_type.num[!is.na(order.tscan)]
  order.tscan.temp = order.tscan[!is.na(order.tscan)]
  cell.names.temp = cell.names[!is.na(order.tscan)]

  kendall.score = cor(cbind(cell_type.num.temp,
                                 order.tscan.temp), method="kendall")[1,2]
  if(kendall.score < 0){
    order.tscan = max(order.tscan, na.rm = T) + 1 - order.tscan
    order.tscan.temp = max(order.tscan.temp) + 1 - order.tscan.temp

    kendall.score = cor(cbind(cell_type.num.temp,
                                   order.tscan.temp), method="kendall")[1,2]
  }
  subpopulation <- data.frame(cell = cell.names.temp,
                              sub = cell_type.num.temp)
  orders = list(order(order.tscan.temp))
  POS.score = orderscore(subpopulation, orders)

  TSCAN.clust = data.clust
  TSCAN.order = order.tscan

  trajectoryTSCAN = plotmclust(TSCAN.clust, show_cell_names = F,
                                    cell_type = cell_type) +
    ggtitle(paste0(name,
                   " (POS: ",round(POS.score,3),
                   ", Cor: ", round(kendall.score,3), ")"))
  pseudoTSCAN = plot_Pseudotime(TSCAN.order, name, cell_type)

  return(list(trajectoryTSCAN = trajectoryTSCAN,
              pseudoTSCAN = pseudoTSCAN,
              POS.score = POS.score,
              kendall.score = kendall.score))
}

plotmclust <- function(mclustobj, x = 1, y = 2, MSTorder = NULL, show_tree = T,
                       show_full_tree = T, show_cell_names = F,
                       cell_name_size = 3, cell_point_size=3, markerexpr = NULL,
                       showcluster = T, cell_type) {
  color_by = "State"
  State = as.factor(paste0("State ", mclustobj$clusterid))

  lib_info_with_pseudo <- data.frame(State=State,sample_name=names(mclustobj$clusterid))
  lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
  S_matrix <- mclustobj$pcareduceres
  pca_space_df <- data.frame(S_matrix[,c(x, y)])
  colnames(pca_space_df) <- c("PC 1","PC 2")
  pca_space_df$`Cell Type` = cell_type
  pca_space_df$sample_name <- row.names(pca_space_df)
  edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
  edge_df$Marker <- markerexpr[edge_df$sample_name]

  g <- ggplot(data = edge_df, aes(x = `PC 1`, y = `PC 2`))

  # `Cell Type` = cell_type
  # g <- g + geom_point(aes(color = `Cell Type`, shape = State), na.rm = TRUE, size = cell_point_size)
  g <- g + geom_point(aes(color = `Cell Type`), size=I(1.5), na.rm = TRUE, size = cell_point_size)

  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }

  if (show_tree) {
    clucenter <- mclustobj$clucenter[,c(x,y)]
    clulines <- NULL
    if (show_full_tree) {
      alledges <- as.data.frame(igraph::get.edgelist(mclustobj$MSTtree),stringsAsFactors=F)
      alledges[,1] <- as.numeric(alledges[,1])
      alledges[,2] <- as.numeric(alledges[,2])
      for (i in 1:nrow(alledges)) {
        clulines <- rbind(clulines, c(clucenter[alledges[i,1],],clucenter[alledges[i,2],]))
      }
    } else {
      if (is.null(MSTorder)) {
        clutable <- table(State)
        alldeg <- igraph::degree(mclustobj$MSTtree)
        allcomb <- expand.grid(as.numeric(names(alldeg)[alldeg == 1]), as.numeric(names(alldeg)[alldeg == 1]))
        allcomb <- allcomb[allcomb[, 1] < allcomb[, 2], ]
        numres <- t(apply(allcomb, 1, function(i) {
          tmp <- as.vector(igraph::get.shortest.paths(mclustobj$MSTtree,
                                                      i[1], i[2])$vpath[[1]])
          c(length(tmp), sum(clutable[tmp]))
        }))
        optcomb <- allcomb[order(numres[, 1], numres[, 2], decreasing = T)[1], ]
        MSTorder <- igraph::get.shortest.paths(mclustobj$MSTtree, optcomb[1],
                                               optcomb[2])$vpath[[1]]
      }
      for (i in 1:(length(MSTorder)-1)) {
        clulines <- rbind(clulines, c(clucenter[MSTorder[i],],clucenter[MSTorder[i+1],]))
      }
    }
    clulines <- data.frame(x=clulines[,1],xend=clulines[,3],y=clulines[,2],yend=clulines[,4])
    g <- g + geom_segment(aes_string(x="x",xend="xend",y="y",yend="yend",size=NULL),data=clulines,size=1)

    # clucenter <- data.frame(x=clucenter[,1],y=clucenter[,2],id=1:nrow(clucenter))
    # g <- g + geom_text(aes_string(label="id",x="x",y="y",size=NULL),data=clucenter,size=10)

  }
  g <- g + guides(colour = guide_legend(override.aes = list(size=5))) +
    xlab(paste0("PC ",x)) + ylab(paste0("PC ",y)) +
    scale_color_tableau() +
    # ggtitle(name) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 14),
          legend.title=element_blank())
  return(g)
}

plot_Pseudotime <- function(order.tscan, name, cell_type){
  df = data.frame(order.tscan, cell_type)

  pTSCAN = ggplot(df, aes(x = order.tscan, y = cell_type, colour = cell_type)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() +
    ggtitle(name) +
    theme_classic() +
    xlab("TSCAN pseudotime") + ylab("Timepoint") +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 14),
          legend.title=element_blank())

  return(pTSCAN)
}
