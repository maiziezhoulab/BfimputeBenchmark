#------------------------------------------------------------------------------#
sra.run.info = read.csv("./SraRunInfo.csv", header = T)
sra.result = read.csv("./sra_result.csv", header = T)
sra.result.extitle = sra.result$Experiment.Title

sra.cell.names = sapply(seq(sra.result.extitle),function(n){
  temp = unlist(strsplit(as.character(sra.result.extitle[n]), " "))
  temp.name = temp[grep("[[]",temp)]
  temp.name = gsub("[[]", "", temp.name)
  temp.name = gsub("[]]", "", temp.name)
  temp.name = gsub("[;]", "", temp.name)
  temp.name = gsub("[.]", "_", temp.name)
  return(temp.name)})

sra.result = data.frame(sra.result, cell.name = sra.cell.names)
sra.inter = intersect(sra.result$Experiment.Accession, sra.run.info$Experiment)
sra.result.match = match(sra.inter, sra.result$Experiment.Accession)
sra.run.info.match = match(sra.inter, sra.run.info$Experiment)

sra.result = sra.result[sra.result.match,]
sra.run.info = sra.run.info[sra.run.info.match, ]

index.df = data.frame(run.name = sra.run.info$Run, cell.name = sra.result$cell.name)
write.csv(index.df, "./indexdf.csv")


#------------------------------------------------------------------------------#
reference = "GRCh37"
outdir = "../count_matrix/"
dir.create(file.path(out_dir), recursive = TRUE)
SRRdir = "./"

index.df = read.csv("./indexdf.csv", row.names = 1, header = T)

create_matrix = function(index_file, count_dir, reference, index.df, outdir, outname){
  KTN = as.character(as.matrix(read.table(index_file)))
  ematrix_raw = sapply(KTN, function(n){
    temp = read.table(paste0(count_dir,"/",n,"_",reference,".count"))
    temp = temp[1:(nrow(temp)-5), 2]
    return(temp)
  })
  genes = read.table(paste0(count_dir,"/", KTN[1],"_",reference,".count"))
  genes = genes[1:(nrow(genes)-5), 1]
  rownames(ematrix_raw) = genes
  colnames(ematrix_raw) = index.df$cell.name[which(index.df$run.name %in% KTN)]
  saveRDS(ematrix_raw, paste0(outdir,"/ematrix_raw_", outname, ".rds"))
  return(0)
}

run.files = list.files(SRRdir)
for(i in seq(run.files)){
  index_file = paste0(SRRdir, run.files[i])
  if(!grepl(".txt",run.files[i])) next
  name_file = gsub(".txt", "", run.files[i])
  count_dir = paste0("../count_", name_file)
  create_matrix(index_file, count_dir, reference, index.df, outdir, name_file)
}


#------------------------------------------------------------------------------#
#--------------------------------------QC--------------------------------------#
#------------------------------------------------------------------------------#
files = list.files(outdir)

file.list.sc = setdiff(files[grep("sc", files)], c("sc_raw.rds", "sc_raw_symbol.rds"))
sc = c()
for(i in 1 : 6) sc = cbind(sc, readRDS(file.list.sc[i]))
colnames(sc) = gsub("00hb4s", "00h", colnames(sc))
sc = sc[, order(colnames(sc))]
saveRDS(sc, paste0(outdir,"/sc_raw.rds"))

file.list.bulk = setdiff(files[grep("bulk", files)], c("bulk_raw.rds", "bulk_raw_symbol.rds"))
bulk = c()
for(i in 1:5) bulk = cbind(bulk, readRDS(file.list.bulk[i]))
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
library(sva)

dir.create(file.path("time_new"), recursive = TRUE)

# SI
sc.SI = sc[, setdiff(seq(colnames(sc)), grep("00h", colnames(sc)))]
sc.SI = sc.SI[rowSums(sc.SI)!=0, ]
sc.SI = preprocess(sc.SI)

bulk.SI = bulk[rowSums(bulk)!=0, ]
bulk.SI = preprocess(bulk.SI)

gene.inter = intersect(rownames(sc.SI), rownames(bulk.SI))
sc.SI = sc.SI[gene.inter, ]
bulk.SI = bulk.SI[gene.inter, ]

bind.matrix = cbind(sc.SI, bulk.SI)
batch = c(rep(1,ncol(sc.SI)), rep(2,ncol(bulk.SI)))
adjusted.matrix = ComBat_seq(bind.matrix, batch=batch, group=NULL)

sc.SI = adjusted.matrix[,seq(ncol(sc.SI))]
bulk.SI = adjusted.matrix[,-(seq(ncol(sc.SI)))]

# write.csv(sc.SI, "./time_new/time_new_SI_qc.csv")
# write.csv(bulk.SI, "./time_new/time_new_bulk_SI_qc.csv")
dir.create(file.path("../time_new/"), recursive = TRUE)
saveRDS(sc.SI, "../time_new/time_new_SI_qc.rds")
saveRDS(bulk.SI, "../time_new/time_new_bulk_SI_qc.rds")

cell_type_SI = sapply(seq(colnames(sc.SI)),function(n){
  temp = unlist(strsplit(colnames(sc.SI)[n],"_"))[2]})

# write.csv(cell_type_SI, "./time_new/time_new_cell_type_SI.csv")
saveRDS(cell_type_SI, "../time_new/time_new_cell_type_SI.rds")





























