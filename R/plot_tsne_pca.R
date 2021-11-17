#' Scatter plot using t-SNE and first two PCs from PCA
#'
#' @param ematrix Expression matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param name Title of the plot.
#' @param cell_type Cell_type of the rows of \code{ematrix}, aka the colour with which
#' you want to paint the dots.
#'
#' @return a list containing PCA plot and t-SNE plot (both are ggplot2 objects)
#' @export
#' @import SingleCellExperiment scater ggplot2 cowplot gridExtra ggthemes
#' @author Zi-Hang Wen \email{wenzihang0506@gmail.com}
#'
#' @examples
#'
plot_tsne_pca <- function(ematrix, name, cell_type){
  res = ematrix
  colnames(res) <- NULL
  rownames(res) <- NULL
  res <- SingleCellExperiment(
    assays = list(logcounts = log2(as.matrix(res) + 1)),
    colData = data.frame(cell_type = cell_type)
  )
  set.seed(123456)
  res <- runPCA(res)
  set.seed(123456)
  res <- runTSNE(res)

  p = plotPCA(res) +
    geom_point(aes(colour = cell_type),size = 1) +
    scale_color_tableau() +
    ggtitle(name) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 12),
          legend.title=element_blank())

  t = plotTSNE(res) +
    geom_point(aes(colour = cell_type),size = 1) +
    scale_color_tableau() +
    ggtitle(name) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 12),
          legend.title=element_blank())
  # + scale_color_brewer(palette = "Paired")

  return(list(p.pca = p,p.tsne = t))
}
