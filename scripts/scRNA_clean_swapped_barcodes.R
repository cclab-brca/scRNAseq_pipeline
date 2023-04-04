#!/home/bioinformatics/software/R/R-4.2.2-gcc-8.5.0/bin/Rscript
require(data.table, quietly = T)
require(DropletUtils, quietly = T)
require(Matrix, quietly = T)

args = commandArgs(trailingOnly = T)

run_dir = args[1]
min_frac = as.numeric(args[2])

mol_info_fns = list.files(path=run_dir, pattern="molecule_info.h5", full.names=T, recursive=T)
x = swappedDrops(mol_info_fns, get.swapped = T, min.frac=min_frac)

names(x$swapped) = gsub(paste0(run_dir, "/"), "", gsub("/outs/molecule_info.h5", "", mol_info_fns))
names(x$cleaned) = names(x$swapped)

df = NULL
for (sample_name in names(x$swapped)) {
  cleaned = x$cleaned[[sample_name]]
  swapped = x$swapped[[sample_name]]
  colnames(cleaned) = colnames(swapped) = paste0(colnames(cleaned), "-1")
  
  bcs = fread(sprintf("%s/%s/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", run_dir, sample_name), header=F)$V1
  gs = fread(sprintf("%s/%s/outs/filtered_feature_bc_matrix/features.tsv.gz", run_dir, sample_name), header=F)
  stopifnot(all(rownames(cleaned) == gs$V1) && all(bcs %in% colnames(cleaned)))
  
  sh_bcs = intersect(bcs, colnames(cleaned))
  
  orig_dir = sprintf("%s/%s/outs/filtered_feature_bc_matrix", run_dir, sample_name)
  
  cleaned_dir = sprintf("%s_swappedDrops_cleaned/", orig_dir)
  dir.create(cleaned_dir, showWarnings = F)
  Matrix::writeMM(cleaned[, sh_bcs], sprintf("%s/matrix.mtx", cleaned_dir))
  write.table(sh_bcs, file = paste0(cleaned_dir, "/barcodes.tsv"), quote = F, row.names = F, col.names = F) 
  file.copy(from = sprintf("%s/features.tsv.gz", orig_dir), to = sprintf("%s/features.tsv.gz", cleaned_dir), overwrite = T)
  xx = lapply(paste(cleaned_dir, c("matrix.mtx", "barcodes.tsv"), sep="/"), R.utils::gzip, overwrite=T)
  
  swapped_dir = sprintf("%s_swappedDrops_swapped/", orig_dir)
  dir.create(swapped_dir, showWarnings = F)
  #fwrite_mm(swapped[, sh_bcs], fname=sprintf("%s/matrix.mtx", swapped_dir), row.names=paste0(swapped_dir, "/features.tsv"), col.names=paste0(swapped_dir, "/barcodes.tsv"))
  Matrix::writeMM(swapped[, sh_bcs], sprintf("%s/matrix.mtx", swapped_dir))
  write.table(sh_bcs, file = paste0(swapped_dir, "/barcodes.tsv"), quote = F, row.names = F, col.names = F) 
  file.copy(from = sprintf("%s/features.tsv.gz", orig_dir), to = sprintf("%s/features.tsv.gz", swapped_dir), overwrite = T)
  xx = lapply(paste(swapped_dir, c("matrix.mtx", "barcodes.tsv"), sep="/"), R.utils::gzip, overwrite=T)
  
  c_df = data.frame(SampleName = sample_name, tot_cleaned=sum(cleaned), tot_swapped=sum(swapped))
  c_df$f_swapped = c_df$tot_swapped / (c_df$tot_swapped + c_df$tot_cleaned)
  
  df = rbind(df, c_df)
}

write.csv(df, sprintf("%s/report/swappedDrops_counts.csv", run_dir), row.names=F)

file.create(sprintf("%s/cleaned_swapped_barcodes.done", run_dir))
