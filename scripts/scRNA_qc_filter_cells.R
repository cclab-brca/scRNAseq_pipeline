#!/home/bioinformatics/software/R/R-4.0.4-gcc-9.2.0/bin/Rscript

# Load a cellranger count generated matrix into Metacell and create a QC passed barcodes list in 
# the outs folder, filtering by non-Mito UMIs per cell and total mito umis ('auto' for automatic cutoffs based on the data).
# 
# If CellBender was used to remove ambient noise, will load its output matrix and generate some comparison plots.

require(data.table, quietly = T)
require(metacell, quietly = T)
require(dplyr, quietly = T)
require(tibble, quietly = T)
require(Matrix, quietly = T)
require(Seurat, quietly = T)

args = commandArgs(trailingOnly = T)
if (length(args) == 0) {
  message("Usage: scRNA_qc_filter_cells.R <base dir> <sample name> <min hg umis> <max hg %mito> <min mm umis> <max mm %mito> <human gene pref> <mouse gene pref>\n")
  stop()
}

mm_col = 'purple'

base_dir = args[1]
sample_name = args[2]
hg_pref = args[7]
mm_pref = args[8]

mito_filter_by = args[9]

qc_csv_fn = sprintf("%s/report/tmp/%s_qcs.csv", base_dir, sample_name)

stopifnot(hg_pref != "NA" || mm_pref != "NA")
hg_pref = ifelse(mm_pref == "NA", "", hg_pref)
mm_pref = ifelse(hg_pref == "NA", "", mm_pref)


if (!mito_filter_by %in% c('frac', 'frac_and_count')) {
  stop(sprintf("Error: unknown MITO_FILTER_BY value (%s), expecting 'frac' or 'frac_and_count'.", mito_filter_by))
}

samp_dir = sprintf("%s/%s/outs", base_dir, sample_name)

scdb_dir = paste0(base_dir, "/scdb") 
figs_dir = paste0(base_dir, "/figs")

dir.create(scdb_dir, showWarnings = F)
dir.create(figs_dir, showWarnings = F)
scdb_init(scdb_dir)
scfigs_init(figs_dir)

# Load counts matrices ----
message("building mat...")
cb_h5_fn = sprintf("%s/cellbender_output_filtered.h5", samp_dir)
cellbender = file.exists(cb_h5_fn)
cleaned_swapped_idir = sprintf("%s/filtered_feature_bc_matrix_swappedDrops_cleaned", samp_dir)
cleaned_swapped_bcs = file.exists(cleaned_swapped_idir)

mat = scmat_read_scmat_10x(matrix_fn = sprintf("%s/filtered_feature_bc_matrix/matrix.mtx.gz", samp_dir), 
                           genes_fn  = sprintf("%s/filtered_feature_bc_matrix/features.tsv.gz", samp_dir),
                           cells_fn  = sprintf("%s/filtered_feature_bc_matrix/barcodes.tsv.gz", samp_dir), 
                           dataset_id = paste0(sample_name, ifelse(cellbender, "_orig", "")))
mat@cell_metadata$seq_batch_id = basename(base_dir)
mat@cell_metadata$amp_batch_id = sample_name
scdb_add_mat(paste0(sample_name, ifelse(cellbender, "_orig", "")), mat)

if (cellbender && cleaned_swapped_bcs) {
  stop(sprintf("Error: both CellBender and swappedDrops were run on this sample, when only one matrix can be used.\nPlease either delete/rename %s or %s to conntinue with only one of them", cb_h5_fn, cleaned_swapped_idir))  
}

# CellBender was used to remove ambient noise. Load both original and cleaned count matrices. ----
if (cellbender) {
  message("Loading count matrix cleaned by CellBender...")
  mf = Seurat::Read10X_h5(cb_h5_fn, unique.features=F)
  empty_cells = colSums(mf) == 0
  if (sum(empty_cells) > 0) {
    message(sprintf("Removing %d cells with no UMIs after ambient noise cleaning.", sum(empty_cells)))
    mf = mf[, !empty_cells]
  }  
  
  unique_genes = names(which(table(rownames(mf)) == 1))
  non_unique_genes = names(which(table(rownames(mf)) > 1))
  mf2 = as.matrix(mf[rownames(mf) %in% non_unique_genes, ])
  message(sprintf("summing up total of %d paralog genes into %d unique genes",
                  nrow(mf2), length(non_unique_genes)))
  mf2s = apply(mf2, 2, function(x) { tapply(x, INDEX = rownames(mf2), FUN = sum) })
  mf = rbind(mf[unique_genes, ], mf2s)
  mf = mf[mat@genes, ]
  
  mf_md = data.frame(row.names=colnames(mf))
  mf_md$amp_batch_id = mf_md$batch_set_id = sample_name
  mf_md$seq_batch_id = basename(base_dir)
  
  scdb_add_mat(sample_name, tgScMat(mf, "umi", mf_md))
  mat = scdb_mat(sample_name)
  
} 

if (cleaned_swapped_bcs) {
  message("loading the swapped barcodes cleaned matrix...")
  mat = scmat_read_scmat_10x(matrix_fn = sprintf("%s/matrix.mtx.gz", cleaned_swapped_idir), 
                             genes_fn  = sprintf("%s/features.tsv.gz", cleaned_swapped_idir),
                             cells_fn  = sprintf("%s/barcodes.tsv.gz", cleaned_swapped_idir), 
                             dataset_id = sample_name)
  mat@cell_metadata$seq_batch_id = basename(base_dir)
  mat@cell_metadata$amp_batch_id = sample_name
  scdb_add_mat(sample_name, mat)
}

mm_uc = colSums(mat@mat[grep(paste0("^", mm_pref), mat@genes), ])
mm_mt = colSums(mat@mat[grep(sprintf("^%smt-", ifelse(hg_pref == "NA", "", paste0(mm_pref, "_"))), mat@genes), ]) 
mm_non_mt = mm_uc - mm_mt
mm_f_mt = ifelse(mm_uc > 0, mm_mt / mm_uc, 0)
mm_n_genes = colSums(mat@mat[grep(paste0("^", mm_pref), mat@genes), ] > 0)

hg_uc = colSums(mat@mat[grep(paste0("^", hg_pref), mat@genes), ])
hg_mt = colSums(mat@mat[grep(sprintf("^%sMT-", ifelse(mm_pref == "NA", "", paste0(hg_pref, "_"))), mat@genes), ]) 
hg_non_mt = hg_uc - hg_mt
hg_f_mt = ifelse(hg_uc > 0, hg_mt / hg_uc, 0)
hg_n_genes = colSums(mat@mat[grep(paste0("^", hg_pref), mat@genes), ] > 0)

is_hg = hg_uc >= mm_uc

# Get %mito and tot non-mito umis cutoffs 
get_bimodal_mid_nadir_cutoff = function(x, from_q=0.2, to_q=0.99, nbins=20) {
  if (length(x) == 0) {
    return(NA)
  }
  from = floor(quantile(x, from_q) * 20) / 20
  to = ceiling(quantile(x, to_q) * 20) / 20
  h = hist(pmin(pmax(x, from), to), breaks=seq(from, to, len=nbins), plot=F)
  h$mids[which.min(h$counts)]
}

# Find qc cutoffs
max_hg_f_mito = ifelse(args[4] == 'auto', get_bimodal_mid_nadir_cutoff(hg_f_mt[is_hg]), ifelse(args[4] == 'from_table', NA, as.numeric(args[4])))
max_mm_f_mito = ifelse(args[6] == 'auto', get_bimodal_mid_nadir_cutoff(mm_f_mt[!is_hg]), ifelse(args[6] == 'from_table', NA, as.numeric(args[6])))

min_hg_non_mt_umis = round(ifelse(args[3] == 'auto', 2 ^ get_bimodal_mid_nadir_cutoff(log2(hg_non_mt[is_hg]), to_q=0.5, nbins=40), ifelse(args[3] == 'from_table', NA, as.numeric(args[3]))))
min_mm_non_mt_umis = round(ifelse(args[5] == 'auto', 2 ^ get_bimodal_mid_nadir_cutoff(log2(mm_non_mt[!is_hg]), to_q=0.5, nbins=40), ifelse(args[5] == 'from_table', NA, as.numeric(args[5]))))

if (any(args[3:6] == 'from_table') && !file.exists(qc_csv_fn)) {
  stop(sprintf("Error: requested to get some QC cutoffs from sample csv but that file does not exist (%s).", qc_csv_fn))
}

if (any(args[3:6] == 'from_table')) {
  message(sprintf("Using cutoffs from %s table for: %s", qc_csv_fn, paste0(c('hg nonMT umis', 'hg fMT', 'mm nonMT umis', 'mm fMT')[args[3:6] == 'from_table'], collapse=", ")))
  qc_csv = read.csv(qc_csv_fn)
  max_hg_f_mito = ifelse(args[4] == 'from_table', qc_csv$hg_max_f_MT_cutoff, max_hg_f_mito)
  max_mm_f_mito = ifelse(args[6] == 'from_table', qc_csv$mm_max_f_MT_cutoff, max_mm_f_mito)
  min_hg_non_mt_umis = ifelse(args[3] == 'from_table', qc_csv$hg_min_nonMT_umi_cutoff, min_hg_non_mt_umis)
  min_mm_non_mt_umis = ifelse(args[5] == 'from_table', qc_csv$mm_min_nonMT_umi_cutoff, min_mm_non_mt_umis)
}
  
max_hg_mito_umis = round(quantile(hg_mt[is_hg & hg_f_mt >= max_hg_f_mito], 0.01))
max_mm_mito_umis = round(quantile(mm_mt[!is_hg & mm_f_mt >= max_mm_f_mito], 0.01))

message(sprintf("%s cutoffs:\n\t\tfMT\t\tMT umis\t\tnon-MT umis\nHuman\t\t%.2f\t\t%d\t\t%d\nMouse\t\t%.2f\t\t%d\t\t%d\n\n",
                sample_name, max_hg_f_mito, max_hg_mito_umis, min_hg_non_mt_umis, max_mm_f_mito, max_mm_mito_umis, min_mm_non_mt_umis))

message(sprintf("Filtering cells by UMI count and mito UMIS (%s)", mito_filter_by))
hg_mt_valid = switch(mito_filter_by,
                     frac = hg_f_mt < max_hg_f_mito,
                     frac_and_count = hg_mt < max_hg_mito_umis | hg_f_mt < max_hg_f_mito)
hg_umi_valid = hg_non_mt >= min_hg_non_mt_umis

mm_mt_valid = switch(mito_filter_by,
                     frac = mm_f_mt < max_mm_f_mito,
                     frac_and_count = mm_mt < max_mm_mito_umis | mm_f_mt < max_mm_f_mito)
mm_umi_valid = mm_non_mt >= min_mm_non_mt_umis

is_valid = NULL
if (hg_pref != "NA" && mm_pref != "NA") {
  is_valid = switch(mito_filter_by,
                    frac = ifelse(is_hg, hg_umi_valid & hg_mt_valid, mm_umi_valid & mm_mt_valid), 
                    frac_and_count = hg_mt_valid & mm_mt_valid & ifelse(is_hg, hg_umi_valid, mm_umi_valid))
  message(sprintf("%d of %d human (%.2f) and %d of %d mouse (%.2f%%) cells passed qc cutoffs.\n\n", sum(is_hg & is_valid), sum(is_hg), 100*mean(is_valid[is_hg]), sum(!is_hg & is_valid), sum(!is_hg), 100*mean(is_valid[!is_hg])))
} else if (hg_pref == "NA") { 
  is_valid = mm_mt_valid & mm_umi_valid 
  message(sprintf("%d of %d mouse cells (%.2f) passed qc cutoffs.\n\n", sum(is_valid), length(is_valid), mean(is_valid)))
} else {
  is_valid = hg_mt_valid & hg_umi_valid
  message(sprintf("%d of %d human cells (%.2f) passed qc cutoffs.\n\n", sum(is_valid), length(is_valid), mean(is_valid)))
}

qc_csv = data.frame(Sample.name = sample_name,
                    f_human = mean(is_hg))

if (hg_pref != "NA") {
  mat@cell_metadata = cbind(mat@cell_metadata,
                          data.frame(hg_tot_umis = hg_uc,
                                     hg_tot_mt_umis = hg_mt,
                                     hg_tot_non_mt_umis = hg_non_mt,
                                     hg_tot_genes = hg_n_genes,
                                     hg_f_mito = hg_f_mt,
                                     hg_mt_good = hg_mt < max_hg_mito_umis | hg_f_mt < max_hg_f_mito,
                                     hg_non_mt_good = hg_non_mt >= min_hg_non_mt_umis))
  qc_csv = qc_csv %>% 
    mutate(n_valid_hg = sum(is_valid & is_hg),
           f_valid_hg = mean(is_valid[is_hg]),
           valid_hg_med_tot_umis = median(hg_uc[is_valid & is_hg]),
           valid_hg_med_nonMT_umis = median(hg_non_mt[is_valid & is_hg]),
           valid_hg_med_genes = median(hg_n_genes[is_valid & is_hg]),
           hg_min_nonMT_umi_cutoff = min_hg_non_mt_umis,
           hg_max_f_MT_cutoff = max_hg_f_mito,
           hg_f_high_MT = mean(!hg_mt_valid[is_hg]),
           hg_high_MT_med_fMT = median(hg_f_mt[!hg_mt_valid & is_hg]))
}
if (mm_pref != "NA") {
  mat@cell_metadata = cbind(mat@cell_metadata,
                            data.frame(mm_tot_umis = mm_uc,
                                       mm_tot_mt_umis = mm_mt,
                                       mm_tot_non_mt_umis = mm_non_mt,
                                       mm_tot_genes = mm_n_genes,
                                       mm_f_mito = mm_f_mt,
                                       mm_mt_good = mm_mt < max_mm_mito_umis | mm_f_mt < max_mm_f_mito,
                                       mm_non_mt_good = mm_non_mt >= min_mm_non_mt_umis))
  qc_csv = qc_csv %>% 
    mutate(n_valid_mm = sum(is_valid & !is_hg),
           f_valid_mm = mean(is_valid[!is_hg]),
           valid_mm_med_tot_umis = median(mm_uc[is_valid & !is_hg]),
           valid_mm_med_nonMT_umis = median(mm_non_mt[is_valid & !is_hg]),
           valid_mm_med_genes = median(mm_n_genes[is_valid & !is_hg]),
           mm_min_nonMT_umi_cutoff = min_mm_non_mt_umis,
           mm_max_f_MT_cutoff = max_mm_f_mito,
           mm_f_high_MT = mean(!mm_mt_valid[!is_hg]),
           mm_high_MT_med_fMT = median(mm_f_mt[!mm_mt_valid & !is_hg]))
}
mat@cell_metadata$valid_by_umis = is_valid

if (cellbender) {
  omat = scdb_mat(paste0(sample_name, "_orig"))
  stopifnot(all(omat@genes == mat@genes))
  shc = intersect(omat@cells, mat@cells)
  
  dm = omat@mat[, shc] - mat@mat[, shc]
  
  f_amb = colSums(dm) / colSums(omat@mat[, shc])
  f_hg_amb = colSums(dm[grep(paste0("^", ifelse(mm_pref == "NA", "", paste0(hg_pref, "_"))), rownames(dm)), ]) / colSums(dm)

  cb_csv = read.csv(qc_csv_fn)
  qc_csv = qc_csv %>% left_join(cb_csv)
  
  qc_csv = qc_csv %>%
    mutate(hg_med_f_amb_noise = median(f_amb[is_hg[shc]]),
           hg_med_f_hg_of_amb_noise = median(f_hg_amb[is_hg[shc]]),
           mm_med_f_amb_noise = median(f_amb[!is_hg[shc]]),
           mm_med_f_hg_of_amb_noise = median(f_hg_amb[!is_hg[shc]]))
  
  mat@cell_metadata[shc, 'f_amb'] = f_amb
  mat@cell_metadata[shc, 'f_hg_g_of_amb'] = f_hg_amb
  
} else {
  qc_csv = qc_csv %>%
    mutate(CellBender_Ncells=NA,
           CellBender_Ndroplets=NA,
           CellBender_FPR=NA,
           CellBender_Epochs=NA, 
           hg_med_f_amb_noise = NA,
           hg_med_f_hg_of_amb_noise = NA,
           mm_med_f_amb_noise = NA,
           mm_med_f_hg_of_amb_noise = NA)
}

scdb_add_mat(sample_name, mat)

md = mat@cell_metadata

# write filtered barcodes and QC metrics
fwrite(as.data.frame(rownames(md)[md$valid_by_umis]), file=sprintf("%s/qc_filt_barcodes.tsv.gz", samp_dir), quote=F, sep="\t", col.names = F, compress="gzip")
write.csv(qc_csv, qc_csv_fn, quote=F, row.names = F)

# some plots ----
if (hg_pref != "NA" && mm_pref != "NA" && sum(is_hg) > 20 && sum(!is_hg) > 20) {
  png(scfigs_fn(sample_name, "all_UMIs_by_specie"), 400, 400)
  pm = par("mai")
  mmd = density(log2(1 + md$mm_tot_umis[!is_hg]))
  hgd = density(log2(1 + md$hg_tot_umis[is_hg]))
  plot(mmd, xlab="UMIs (log2)", main=sprintf("%s All", sample_name), lwd=2, col=mm_col, xlim=range(c(mmd$x, hgd$x)), ylim=range(c(mmd$y, hgd$y)))
  lines(hgd, lwd=2, col='black')
  legend("topleft", legend=c("Human", "Mouse"), lwd=2, col=c('black', mm_col), bty='n')
  dev.off()

  png(scfigs_fn(sample_name, "nonMT_UMIs_by_specie"), 400, 400)
  mmd = density(log2(1 + md$mm_tot_non_mt_umis[!is_hg]))
  hgd = density(log2(1 + md$hg_tot_non_mt_umis[is_hg]))
  plot(mmd, xlab="UMIs (log2)", main=sprintf("%s non-MT", sample_name), lwd=2, col=mm_col, xlim=range(c(mmd$x, hgd$x)), ylim=range(c(mmd$y, hgd$y)))
  lines(hgd, lwd=2, col='black')
  abline(v=log2(c(min_hg_non_mt_umis, min_mm_non_mt_umis)), lty=2, col=c('black', mm_col))
  legend("topleft", legend=c("Human", "Mouse"), lwd=2, col=c('black', mm_col), bty='n')
  dev.off()

  png(scfigs_fn(sample_name, "fMito_vs_logUMIs"), 1000, 400)
  layout(matrix(1:4, nrow=1), widths=c(4,1,4,1))
  par(cex=1)
  par(mai=c(pm[1], pm[2]+pm[4]-0.1, pm[3], 0.1))
  plot(log2(md$hg_tot_non_mt_umis[is_hg]), md$hg_f_mito[is_hg], pch=19, cex=0.3, col=ifelse(md$valid_by_umis[is_hg], 'black', 'red'), xlab="Non-MT UMIs (log2)", ylab="%mito", main=sprintf("%s: Human (%d)", sample_name, sum(md$valid_by_umis[is_hg])), ylim=0:1)
  abline(v=log2(min_hg_non_mt_umis), lty=2)
  abline(h=max_hg_f_mito, lty=2)
  h = hist(md$hg_f_mito[is_hg], breaks=seq(0, 1, len=50), plot=F)
  par(mai=c(pm[1], 0, pm[3], 0.1))
  barplot(h$count, space=0, col='darkgray', horiz=T, border=NA)
  
  par(mai=c(pm[1], pm[2]+pm[4]-0.1, pm[3], 0.1))
  plot(log2(md$mm_tot_non_mt_umis[!is_hg]), md$mm_f_mito[!is_hg], pch=19, cex=0.3, col=ifelse(md$valid_by_umis[!is_hg], 'black', 'red'), xlab="non-MT UMIs (log2)", ylab="%mito", main=sprintf("%s: Mouse (%d)", sample_name, sum(md$valid_by_umis[!is_hg])), ylim=0:1)
  abline(v=log2(min_mm_non_mt_umis), lty=2)
  abline(h=max_mm_f_mito, lty=2)
  h = hist(md$mm_f_mito[!is_hg], breaks=seq(0, 1, len=50), plot=F)
  par(mai=c(pm[1], 0, pm[3], 0.1))
  barplot(h$count, space=0, col='darkgray', horiz=T, border=NA)
  dev.off()

  png(scfigs_fn(sample_name, "hg_vs_mm_UMIs"), 800, 400)
  par(mfrow=1:2)
  plot(log2(md$hg_tot_umis), log2(md$mm_tot_umis), pch=19, cex=0.3, col=ifelse(md$valid_by_umis, 'black', 'red'), xlab="Human UMIs (log2)", ylab="Mouse UMIs (log2)", main=paste(sample_name, "All"))
  plot(log2(md$hg_tot_non_mt_umis), log2(md$mm_tot_non_mt_umis), pch=19, cex=0.3, col=ifelse(md$valid_by_umis, 'black', 'red'), xlab="Human UMIs (log2)", ylab="Mouse UMIs (log2)", main=paste(sample_name, "non-MT"))
  dev.off()

  if (cellbender) {
    mdcb = md[!is.na(md$f_amb), ]
    cb_is_hg = mdcb$hg_tot_umis > mdcb$mm_tot_umis
    
    png(scfigs_fn(sample_name, "fAmb_by_specie"), 400, 400)
    mmd = density(mdcb$f_amb[!cb_is_hg])
    hgd = density(mdcb$f_amb[cb_is_hg])
    plot(mmd, xlab="% Ambient noise", main=sprintf("%s %%Amb", sample_name), lwd=2, col=mm_col, xlim=range(c(mmd$x, hgd$x)), ylim=range(c(mmd$y, hgd$y)))
    lines(hgd, lwd=2, col='black')
    legend("topright", legend=c("Human", "Mouse"), lwd=2, col=c('black', mm_col), bty='n')
    dev.off()

    png(scfigs_fn(sample_name, "fAmb_vs_umis_by_specie"), 800, 400)
    par(mfrow=1:2)
    plot(log2(mdcb$hg_tot_non_mt_umis[cb_is_hg]), mdcb$f_amb[cb_is_hg], pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis[cb_is_hg], 'black', 'red'), xlab="Human non-MT UMIs (log2)", ylab="% Amb", main=paste(sample_name, "Human"))
    plot(log2(mdcb$mm_tot_non_mt_umis[!cb_is_hg]), mdcb$f_amb[!cb_is_hg], pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis[!cb_is_hg], 'black', 'red'), xlab="Mouse non-MT UMIs (log2)", ylab="% Amb", main=paste(sample_name, "Mouse"))
    dev.off()

    png(scfigs_fn(sample_name, "fMT_vs_Amb_by_specie"), 800, 400)
    par(mfrow=1:2)
    plot(mdcb$f_amb[cb_is_hg], mdcb$hg_f_mito[cb_is_hg], pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis[cb_is_hg], 'black', 'red'), ylab="% Human MT UMIs", xlab="% Amb", main=paste(sample_name, "Human"))
    plot(mdcb$f_amb[!cb_is_hg], mdcb$mm_f_mito[!cb_is_hg], pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis[!cb_is_hg], 'black', 'red'), ylab="% Mouse MT UMIs", xlab="% Amb", main=paste(sample_name, "Mouse"))
    dev.off()
  }

} else if (hg_pref != "NA" && sum(is_hg) > 20 ) {
  png(scfigs_fn(sample_name, "hg_all_UMIs"), 400, 400)
  pm = par("mai")
  plot(density(log2(1 + md$hg_tot_umis)), xlab="UMIs (log2)", main=sprintf("%s All", sample_name), lwd=2)
  legend("topleft", legend="Human", lwd=2, col='black', bty='n')
  dev.off()
  
  png(scfigs_fn(sample_name, "hg_nonMT_UMIs"), 400, 400)
  plot(density(log2(1 + md$hg_tot_non_mt_umis)), xlab="UMIs (log2)", main=sprintf("%s non-MT", sample_name), lwd=2)
  abline(v=log2(min_hg_non_mt_umis), lty=2, col='black')
  legend("topleft", legend="Human", lwd=2, col='black', bty='n')
  dev.off()
  
  png(scfigs_fn(sample_name, "hg_fMito_vs_logUMIs"), 500, 400)
  layout(matrix(1:2, nrow=1), widths=c(4,1))
  par(cex=1)
  par(mai=c(pm[1], pm[2]+pm[4]-0.1, pm[3], 0.1))
  plot(log2(md$hg_tot_umis), md$hg_f_mito, pch=19, cex=0.5, col=ifelse(md$hg_f_mito >= max_hg_f_mito, 'red', 'black'), xlab="UMIs (log2)", ylab="%mito", main=sprintf("%s: Human (%d)", sample_name, nrow(md)), ylim=0:1)
  h = hist(md$hg_f_mito, breaks=seq(0, 1, len=50), plot=F)
  par(mai=c(pm[1], 0, pm[3], 0.1))
  barplot(h$count, space=0, col='darkgray', horiz=T, border=NA)
  dev.off()
  
  if (cellbender) {
    mdcb = md[!is.na(md$f_amb), ]
    
    png(scfigs_fn(sample_name, "hg_fAmb"), 400, 400)
    plot(density(mdcb$f_amb), xlab="% Ambient noise", main=sprintf("%s %%Amb", sample_name), lwd=2)
    legend("topright", legend="Human", lwd=2, col='black', bty='n')
    dev.off()
    
    png(scfigs_fn(sample_name, "hg_fAmb_vs_umis"), 400, 400)
    plot(log2(mdcb$hg_tot_non_mt_umis), mdcb$f_amb, pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis, 'black', 'red'), xlab="Human non-MT UMIs (log2)", ylab="% Amb", main=sample_name)
    dev.off()
    
    png(scfigs_fn(sample_name, "hg_fMT_vs_Amb"), 400, 400)
    plot(mdcb$f_amb, mdcb$hg_f_mito, pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis, 'black', 'red'), ylab="% Human MT UMIs", xlab="% Amb", main=sample_name)
    dev.off()
  }
  
} else if (sum(!is_hg) > 20) {
  png(scfigs_fn(sample_name, "mm_all_UMIs"), 400, 400)
  pm = par("mai")
  plot(density(log2(1 + md$mm_tot_umis)), xlab="UMIs (log2)", main=sprintf("%s All", sample_name), lwd=2, col='red')
  legend("topleft", legend="Mouse", lwd=2, col='red', bty='n')
  dev.off()
  
  png(scfigs_fn(sample_name, "mm_nonMT_UMIs"), 400, 400)
  plot(density(log2(1 + md$mm_tot_non_mt_umis)), xlab="UMIs (log2)", main=sprintf("%s non-MT", sample_name), lwd=2, col='red')
  abline(v=log2(min_mm_non_mt_umis), lty=2, col='red')
  legend("topleft", legend="Mouse", lwd=2, col='red', bty='n')
  dev.off()
  
  png(scfigs_fn(sample_name, "mm_fMito_vs_logUMIs"), 500, 400)
  layout(matrix(1:2, nrow=1), widths=c(4,1))
  par(cex=1)
  par(mai=c(pm[1], pm[2]+pm[4]-0.1, pm[3], 0.1))
  plot(log2(md$mm_tot_umis), md$mm_f_mito, pch=19, cex=0.5, col=ifelse(md$mm_f_mito >= max_mm_f_mito, 'red', 'black'), xlab="UMIs (log2)", ylab="%mito", main=sprintf("%s: Mouse (%d)", sample_name, nrow(md)), ylim=0:1)
  h = hist(md$mm_f_mito, breaks=seq(0, 1, len=50), plot=F)
  par(mai=c(pm[1], 0, pm[3], 0.1))
  barplot(h$count, space=0, col='darkgray', horiz=T, border=NA)
  dev.off()
  
  if (cellbender) {
    mdcb = md[!is.na(md$f_amb), ]
    
    png(scfigs_fn(sample_name, "mm_fAmb"), 400, 400)
    plot(density(mdcb$f_amb), xlab="% Ambient noise", main=sprintf("%s %%Amb", sample_name), lwd=2)
    legend("topright", legend="Mouse", lwd=2, col='black', bty='n')
    dev.off()
    
    png(scfigs_fn(sample_name, "mm_fAmb_vs_umis"), 400, 400)
    plot(log2(mdcb$mm_tot_non_mt_umis), mdcb$f_amb, pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis, 'black', 'red'), xlab="Mouse non-MT UMIs (log2)", ylab="% Amb", main=sample_name)
    dev.off()
    
    png(scfigs_fn(sample_name, "mm_fMT_vs_Amb"), 400, 400)
    plot(mdcb$f_amb, mdcb$mm_f_mito, pch=19, cex=0.3, col=ifelse(mdcb$valid_by_umis, 'black', 'red'), ylab="% Mouse MT UMIs", xlab="% Amb", main=sample_name)
    dev.off()
  }
  
} 
  

