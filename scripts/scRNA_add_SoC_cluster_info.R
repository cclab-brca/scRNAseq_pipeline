#!/home/bioinformatics/software/R/R-4.0.4-gcc-9.2.0/bin/Rscript

# Find optimal K out of the souporcell clustering solutions and add the selected clustering info to the mat

require(Matrix, quietly = T)
require(metacell, quietly = T)
require(dplyr, quietly = T)
require(tgstat, quietly = T)

args = commandArgs(trailingOnly = T)

base_dir = args[1]
sample_name = args[2]
common_variants = args[3]
max_k = as.numeric(args[4])
ll_improve_cutoff = as.numeric(args[5])

idir = sprintf("%s/%s/outs/souporcell_hg", base_dir, sample_name)

scdb_dir = paste0(base_dir, "/scdb") 
figs_dir = paste0(base_dir, "/figs")

scdb_init(scdb_dir)
scfigs_init(figs_dir)

mat = scdb_mat(paste0(sample_name, "_hg"))
stopifnot(!is.null(mat))

# Find best K (first K that K+1 improved log loss by less than ll_improve_cutoff)
mean_LL = sapply(1:max_k, function(k) { 
  x = read.table(sprintf("%s/k%d_comVar%s/clusters.tsv", idir, k, common_variants), header=T) %>% filter(status == 'singlet'); 
  mean(apply(x[,6:ncol(x), drop=F], 1, max)) 
  })

best_k = min(which(mean_LL[-1] / mean_LL[-max_k] >= (1 - ll_improve_cutoff)))

# elbow plot
png(scfigs_fn(sample_name, sprintf("soc_on_hg_cells_comVar%s_maxK%d_fLLcutoff%.2f_LL_vs_K", common_variants, max_k, ll_improve_cutoff)), 300, 300)
plot(1:max_k, mean_LL, pch=ifelse(1:max_k %in% best_k, 19, 21), main=sample_name, xlab='K', ylab="singlets mean log loss")
dev.off()

# add clusters assignments to mat
best_dir = sprintf("%s/k%d_comVar%s", idir, best_k, common_variants)
mat@cell_metadata$soc_cl = NA
cl = read.table(sprintf("%s/clusters.tsv", best_dir), header=T) %>% filter(status == 'singlet')
mat@cell_metadata[cl$barcode, 'soc_cl'] = as.numeric(cl$assignment) + 1
scdb_add_mat(paste0(sample_name, "_hg"), mat)

# compute VAF over SNVs per cluster
vcf = read.table(sprintf("%s/souporcell_merged_sorted_vcf.vcf.gz", best_dir), header=F, sep="\t")

ref = readMM(sprintf("%s/ref.mtx", best_dir))
alt = readMM(sprintf("%s/alt.mtx", best_dir))

bcs = read.table(sprintf("%s/barcodes.tsv", best_dir))$V1

colnames(ref) = colnames(alt) = bcs
rownames(ref) = rownames(alt) = paste(vcf$V1, vcf$V2, vcf$V4, vcf$V5, sep="_")

f = rowSums(ref + alt) > 0

ref = ref[f, ]
alt = alt[f, ]

alt_c = tgs_matrix_tapply(as.matrix(alt[, cl$barcode]), as.numeric(cl$assignment) + 1, sum)
ref_c = tgs_matrix_tapply(as.matrix(ref[, cl$barcode]), as.numeric(cl$assignment) + 1, sum)
vaf_c = alt_c / (alt_c + ref_c)
write.csv(t(vaf_c), sprintf("%s/clusters_vaf.csv", idir))

file.create(sprintf("%s/%s/outs/soc_on_hg_singlets.done", base_dir, sample_name))
