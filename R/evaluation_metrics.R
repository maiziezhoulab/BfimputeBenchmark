#' Four types of evaluation metrics
#'
#' @param ematrix Expression matrix with rows corresponding to genes and
#' columns corresponding to cells.
#' @param true_label True cell labels
#' @param npc Number of principal components used for dimension reduction
#'
#' @return a list containing adjusted Rand Index, Jaccard Index,
#' normalized mutual information, purity score, and silhouette
#' @export
#' @import Spectrum stats cluster
#' @author Zi-Hang Wen \email{wenzihang0506@gmail.com}
#'
#' @examples
#'
evaluation_metrics <- function(ematrix, true_label, npc = 2){
  old_clust = true_label
  Kcluster = sum(!is.na(unique(old_clust)))
  set.seed(Kcluster)
  pca = prcomp(t(log10(ematrix+1)))
  pc2 = t(pca$x[,1:npc])
  spec_res = Spectrum(pc2,method=3,maxk=10,fixk=Kcluster)
  new_clust = spec_res$assignments

  Knumbers = rep(0, Kcluster)
  for(ii in 1:Kcluster){
    Knumbers[ii] = sum(new_clust == ii)
  }
  print("cluster sizes:")
  print(Knumbers)

  nold = sum(!is.na(unique(old_clust)))
  nnew = sum(!is.na(unique(new_clust)))
  nall = sum(!is.na(old_clust))

  on_matrix = matrix(rep(0,nold*nnew),nrow = nold)

  for(i in 1:nall){
    on_matrix[old_clust[i],new_clust[i]] = on_matrix[old_clust[i],new_clust[i]]+1
  }

  a = sum(on_matrix*(on_matrix-1)/2)
  b = sum(apply(on_matrix,2,function(x) 0.5*(sum(x %o% x) - sum(x^2))))
  c = sum(apply(on_matrix,1,function(x) 0.5*(sum(x %o% x) - sum(x^2))))
  d = sum(sapply(1:dim(on_matrix)[1],function(x) sum(sapply(1:dim(on_matrix)[2], function(y) on_matrix[x,y]*on_matrix[-x,-y]))))/2

  #####-------------------------adjusted Rand index-------------------------#####
  ari = (a - (a+b)*(a+c)/(a+b+c+d)) / (0.5*(a+b+a+c) - (a+b)*(a+c)/(a+b+c+d))

  #ari_p = rep(0,3)
  #ari_p[1] = sum(apply(on_matrix,c(1,2),function(x) return(x*(x-1)/2)))
  #ari_p[2] = sum(sapply(rowSums(on_matrix),function(x) return(x*(x-1)/2))) *
  #  sum(sapply(colSums(on_matrix),function(x) return(x*(x-1)/2))) /
  #  (sum(on_matrix)*(sum(on_matrix)-1)/2)
  #ari_p[3] = 1/2 * (sum(sapply(rowSums(on_matrix),function(x) return(x*(x-1)/2))) +
  #                 sum(sapply(colSums(on_matrix),function(x) return(x*(x-1)/2))))
  #ari = (ari_p[1] - ari_p[2])/(ari_p[3] - ari_p[2])

  #####-------------------------Jaccard index-------------------------#####
  ji = a/(a+b+c)

  #####-------------------------normalized mutual information-------------------------#####
  Hold = sum(-rowSums(on_matrix)/sum(on_matrix)*log2(rowSums(on_matrix)/sum(on_matrix)))
  Hnew = sum(-colSums(on_matrix)/sum(on_matrix)*log2(colSums(on_matrix)/sum(on_matrix)))
  Hon = sum(-colSums(on_matrix)/sum(on_matrix)*colSums(
    t(t(on_matrix)/colSums(on_matrix))*log2(t(t(on_matrix)/colSums(on_matrix))), na.rm = TRUE
  ))
  nmi = 2*(Hold-Hon)/(Hold+Hnew)
  #Hon = 0
  #for(i in 1:ncol(on_matrix)){
  #  Hon = Hon -
  #    colSums(on_matrix)[i]/sum(on_matrix) *
  #    sum(on_matrix[,i]/colSums(on_matrix)[i] * log2(on_matrix[,i]/colSums(on_matrix)[i]),na.rm = TRUE)
  #}

  #####-------------------------purity score-------------------------#####
  ps = sum(apply(on_matrix,2,max))/sum(on_matrix)

  #####--------------------------silhouette--------------------------#####
  counts.dist = dist(t(pc2))
  counts.sil = silhouette(new_clust, counts.dist)
  counts.sil.summary = summary(counts.sil)
  sil.all = list(widths = counts.sil[,3], summary = counts.sil.summary)
  return(list(ari = ari, ji = ji, nmi = nmi, ps = ps, silhouette = sil.all))
  #return(c(ari,ji,nmi,ps))
}

