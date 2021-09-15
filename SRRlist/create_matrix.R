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











