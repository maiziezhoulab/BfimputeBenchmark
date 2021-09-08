#' Title Calculate Pearson Correlation between cells from the same cell type and
#' from different cell type. Also calculate the standard deviation
#'
#' @param res Expression count matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param cell_type A vector of cell types with each element representing a cell
#' type for the cell in the corresponding position in \code{res}
#'
#' @return a list containing two vectors: Pearson scores and their std (the
#' order is same cell type and diff cell type)
#' @export
#' @import stats
#' @author Zi-Hang Wen \email{wenzihang0506@gmail.com}
#'
#' @examples
PearsonStd <- function(res, cell_type){
  pearson_all = cor(res)

  pearson = c()
  pearson[1] = mean(sapply(cell_index.C, function(j){
    mean(pearson_all[cell_type==j, cell_type==j])
  }))
  pearson[2] = mean(sapply(cell_index.C, function(j){
    mean(pearson_all[cell_type==j, cell_type!=j])
  }))
  std = c()
  std[1] = mean(sapply(cell_index.C, function(j){
    sd(pearson_all[cell_type==j, cell_type==j])
  }))
  std[2] = mean(sapply(cell_index.C, function(j){
    sd(pearson_all[cell_type==j, cell_type!=j])
  }))
  return(list(pearson = pearson, std = std))
}
