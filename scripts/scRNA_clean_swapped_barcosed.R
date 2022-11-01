#!/home/bioinformatics/software/R/R-4.0.4-gcc-9.2.0/bin/Rscript
require(data.table, quietly = T)
require(dropletUtils, quietly = T)

args = commandArgs(trailingOnly = T)

run_dir = args[1]
min_frac = as.numeric(args[2])

mol_info_fns = list.files(path=run_dir, pattern="molecule_info.h5", full.names=T, recursive=T)
x = swappedDrops(mol_info_fns, get.swapped = T, min.frac=min_frac)

names(x$swapped) = gsub(paste0(run_idir, "/"), "", gsub("/outs/molecule_info.h5", "", fns))
names(x$cleaned) = names(x$swapped)

df = NULL
for (sample_name in names(x$swapped)) {
  cleaned = x$cleaned[[sample_name]]
  swapped = x$swapped[[sample_name]]
  colnames(cleaned) = colnames(swapped) = paste0(colnamed(cleaned), "-1")
  
  bcs = fread(sprintf("%s/%s/outs/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", run_dir, sample_name))
  gs = fread(sprintf("%s/%s/outs/outs/filtered_feature_bc_matrix/features.tsv.gz", run_dir, sample_name), header=F)
  stopifnot(all(rownames(cleaned) == gs) && all(bcs %in% colnames(cleaned)))
  
  sh_bcs = intersect(bcs, colnames(cleaned))
  
  orig_dir = sprintf("%s/%s/outs/outs/filtered_feature_bc_matrix/", run_dir, sample_name)
  
  cleaned_dir = sprintf("%s_swappedDrops_cleaned/", orig_dir)
  dir.create(cleaned_dir, showWarnings = F)
  writeMM(cleaned[, sh_bcs], sprintf("%s/matrix.mtx.gz", cleaned_dir))
  file.symlink(from=paste0(orig_dir, "/features.tsv.gz"), to = paste0(cleaned_dir, "/features.tsv.gz"))
  file.symlink(from=paste0(orig_dir, "/barcodes.tsv.gz"), to = paste0(cleaned_dir, "/barcodes.tsv.gz"))
  
  swapped_dir = sprintf("%s_swappedDrops_swapped/", orig_dir)
  dir.create(swapped_dir, showWarnings = F)
  writeMM(swapped[, sh_bcs], sprintf("%s/matrix.mtx.gz", swapped_dir))
  file.symlink(from=paste0(orig_dir, "/features.tsv.gz"), to = paste0(swapped_dir, "/features.tsv.gz"))
  file.symlink(from=paste0(orig_dir, "/barcodes.tsv.gz"), to = paste0(swapped_dir, "/barcodes.tsv.gz"))
  
  c_df = data.frame(SampleName = sample_name, tot_cleaned=sum(cleaned), tot_swapped=sum(swapped))
  c_df$f_swapped = c_df$tot_swapped / (c_df$tot_swapped + c_df$tot_cleaned)
  
  df = rbind(df, c_df)
}

write.csv(df, sprintf("%s/report/swappedDrops_counts.csv", run_dir), row.names=F)

file.create(sprintf("%s/cleaned_swapped_barcodes.done", run_dir))
