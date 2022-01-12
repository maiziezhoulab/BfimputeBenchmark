#------------------------------------------------------------------------------#
#--------------------------------------QC--------------------------------------#
#------------------------------------------------------------------------------#
outdir = "../count_matrix/"
files = list.files(outdir)

file.list.sc = setdiff(files[grep("sc", files)], c("sc_raw.rds", "sc_raw_symbol.rds"))
sc = c()
for(i in 1 : 6) sc = cbind(sc, readRDS(paste0(outdir, file.list.sc[i])))
colnames(sc) = gsub("00hb4s", "00h", colnames(sc))
sc = sc[, order(colnames(sc))]
saveRDS(sc, paste0(outdir,"/sc_raw.rds"))

file.list.bulk = setdiff(files[grep("bulk", files)], c("bulk_raw.rds", "bulk_raw_symbol.rds"))
bulk = c()
for(i in 1:5) bulk = cbind(bulk, readRDS(paste0(outdir, file.list.bulk[i])))
bulk = bulk[, order(colnames(bulk))]
saveRDS(bulk, paste0(outdir,"/bulk_raw.rds"))


#------------------------------------------------------------------------------#
library(org.Hs.eg.db)
keys.sc = rownames(sc)
list.gene = select(org.Hs.eg.db, keys=keys.sc, columns = c("ENTREZID","SYMBOL"),
                   keytype="ENSEMBL")
list.gene.id = c()
for(i in seq(keys.sc)){
  list.gene.id = c(list.gene.id, which(list.gene == keys.sc[i])[1])
}
list.gene = list.gene[list.gene.id,]

sc = sc[which(!is.na(list.gene$SYMBOL)),]
rownames(sc) = list.gene$SYMBOL[which(!is.na(list.gene$SYMBOL))]
bulk = bulk[which(!is.na(list.gene$SYMBOL)),]
rownames(bulk) = list.gene$SYMBOL[which(!is.na(list.gene$SYMBOL))]

row.order = order(rownames(sc))
sc = sc[row.order, ]
bulk = bulk[row.order, ]
saveRDS(sc, paste0(outdir,"/sc_raw_symbol.rds"))
saveRDS(bulk, paste0(outdir,"/bulk_raw_symbol.rds"))


#------------------------------------------------------------------------------#
# From https://github.com/gongx030/scDatasets/blob/master/R/scDatasets.r
preprocess <- function(X, min.expressed.gene = 2000, min.expressed.cell = 2, max.expressed.ratio = 1, normalize.by.size.effect = FALSE){

  M0 <- ncol(X)
  N0 <- nrow(X)

  cat(sprintf('[%s] number of input genes(nrow(X))=%.d\n', Sys.time(), N0))
  cat(sprintf('[%s] number of input cells(ncol(X))=%.d\n', Sys.time(), M0))
  X <- X[, Matrix::colSums(X > 1) >= min.expressed.gene, drop = FALSE]
  M <- ncol(X)
  cat(sprintf('[%s] number of input cells that express at least %d genes=%.d\n', Sys.time(), min.expressed.gene, M))
  X <- X[Matrix::rowSums(X > 1) <= max.expressed.ratio * M & Matrix::rowSums(X > 1) >= min.expressed.cell, , drop = FALSE]
  N <- nrow(X)
  cat(sprintf('[%s] number of input genes that are expressed in at least %d cells and at most %.0f%% cells=%.d\n', Sys.time(), min.expressed.cell, max.expressed.ratio * 100, N))
  cat(sprintf('[%s] sparsity of expression matrix=%.1f%%\n', Sys.time(), 100 * (N * M - sum(X > 0)) / (N * M)))

  if (normalize.by.size.effect){
    cat(sprintf('[%s] scaling raw read counts by size factor\n', Sys.time()))
    sf <- apply((X + 1) / exp(Matrix::rowMeans(log(X + 1))), 2, median)
    X <- t(t(X) / sf)
  }

  as.matrix(X)

}


#------------------------------------------------------------------------------#
# SI
sc.SI = sc[, setdiff(seq(colnames(sc)), grep("00h", colnames(sc)))]
sc.SI = sc.SI[rowSums(sc.SI)!=0, ]
sc.SI = preprocess(sc.SI)

bulk.SI = bulk[rowSums(bulk)!=0, ]
bulk.SI = preprocess(bulk.SI)

gene.inter = intersect(rownames(sc.SI), rownames(bulk.SI))
sc.SI = sc.SI[gene.inter, ]
bulk.SI = bulk.SI[gene.inter, ]


dir.create(file.path("../time_new/"), recursive = TRUE)
# write.csv(sc.SI, "./time_new/time_new_SI_qc.csv")
# write.csv(bulk.SI, "./time_new/time_new_bulk_SI_qc.csv")
saveRDS(sc.SI, "../time_new/time_new_SI_qc.rds")
saveRDS(bulk.SI, "../time_new/time_new_bulk_SI_qc.rds")

cell_type_SI = sapply(seq(colnames(sc.SI)),function(n){
  temp = unlist(strsplit(colnames(sc.SI)[n],"_"))[2]})

# write.csv(cell_type_SI, "./time_new/time_new_cell_type_SI.csv")
saveRDS(cell_type_SI, "../time_new/time_new_cell_type_SI.rds")








