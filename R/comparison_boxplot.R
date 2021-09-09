#' Box plot of gene expression comparison between single-cell and bulk
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
#' @param gene The marker gene going to be compared between single-cell and bulk
#'
#' @return the ggplot2 boxplot object
#' @export
#' @import ggplot2 ggpubr gridExtra cowplot
#'
#' @examples
#'
comparison_boxplot <- function(bulk, bulk_cell_type,
                               ematrix, sc_cell_type, name,
                               gene){
  counts.bulk = as.matrix(bulk)
  counts.sc = as.matrix(ematrix)

  gene.inter = intersect(rownames(counts.bulk), rownames(counts.sc))
  counts.bulk = counts.bulk[gene.inter, ]
  counts.sc = counts.sc[gene.inter, ]

  colnames(counts.bulk) = paste0("bulk_",colnames(counts.bulk))
  colnames(counts.sc) = gsub("H9.", "H9_", colnames(counts.sc))
  colnames(counts.sc) = paste0("single-cell_",colnames(counts.sc))

  bulk.ct.unique = unique(bulk_cell_type)

  sc.ct.unique = unique(sc_cell_type)
  sc.ct.num = sapply(seq(sc.ct.unique),function(ii){
    length(which(sc_cell_type == sc.ct.unique[ii]))
  })
  sc.ct.num.cum = c(cumsum(c(1, sc.ct.num)))

  bulk_y_temp = sapply(1:5,function(ii){
    temp = list(mean(log10(counts.bulk[gene, bulk_cell_type == bulk.ct.unique[ii]]+1)))
    names(temp) = colnames(counts.bulk[,bulk_cell_type == bulk.ct.unique[ii]])[1]
    return(temp)
  })
  bulk_y_temp = t(as.matrix(unlist(bulk_y_temp)))

  bulk_y_temp.ct = sapply(1:ncol(bulk_y_temp),function(ii){
    unlist(strsplit(colnames(bulk_y_temp)[ii],"_"))[3]})

  bulk_y = rep(NA,length(sc_cell_type))
  for(ii in seq(bulk_y_temp.ct)){
    if(bulk_y_temp.ct[ii]%in%sc.ct.unique){
      bulk_y[sc.ct.num.cum[which(sc.ct.unique == bulk_y_temp.ct[ii])]] = bulk_y_temp[ii]
    }
  }

  Expression = log10(as.matrix(counts.sc[gene,])+1)
  boxdf = data.frame(Expression, cell_type = sc_cell_type, bulk_y)
  picture_list = ggboxplot(boxdf, x = "cell_type",y = "Expression",
                           color = "cell_type",
                           xlab = "Time course",
                           ylab = "log10(Expression+1)",
                           add = "point") +
    xlab(NULL) +
    labs(title = name) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
    scale_y_continuous(limits = c(0, 3.5)) +
    geom_point(aes(x = cell_type, y = bulk_y),shape = 17, size = 4, na.rm = TRUE)
  return(picture_list)
}

