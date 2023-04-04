#!/home/bioinformatics/software/R/R-4.2.2-gcc-8.5.0/bin/Rscript
require(dplyr, quietly=T)

# Collate QC metrics of the run's sample, add them to a master table (base_dir/scRNAseq_samples_metrics.csv) and link cellranger html reports in a single directory.
args = commandArgs(trailingOnly = T)

base_dir = args[1]
run_name = args[2]

params_fn = sprintf("%s/runs/%s/scRNA_params.sh", base_dir, run_name)
stopifnot(file.exists(params_fn))
chem = gsub("CHEM=", "", grep("CHEM", readLines(params_fn), v=T))
stopifnot(chem %in% c('SC3Pv3', 'SC3Pv2'))

samp_md_fn = list.files(sprintf("%s/runs/%s/metadata", base_dir, run_name), pattern = "*csv", full.names = T)
stopifnot(length(samp_md_fn) == 1)
samp_md = read.csv(samp_md_fn, header=T)


qc_metrics = NULL

selected_cr_columns = c("Estimated Number of Cells", "Fraction Reads in Cells", "Mean Reads per Cell", 
                        "Number of Reads", "Q30 Bases in Barcode", "Q30 Bases in RNA Read", "Q30 Bases in UMI", 
                        "Reads Mapped Antisense to Gene", "Reads Mapped Confidently to Exonic Regions", 
                        "Reads Mapped Confidently to Genome", "Reads Mapped Confidently to Intergenic Regions", 
                        "Reads Mapped Confidently to Intronic Regions", "Reads Mapped Confidently to Transcriptome", 
                        "Reads Mapped to Genome", "Sequencing Saturation", "Valid Barcodes")

for (sample_name in samp_md$Sample.name) {
  message(sample_name)
  
  # cellranger htmls
  to_fn = sprintf("%s/runs/%s/report/%s_cellranger_summary.html", base_dir, run_name, sample_name)
  if (!file.exists(to_fn)) {
    file.copy(sprintf("%s/runs/%s/%s/outs/web_summary.html", base_dir, run_name, sample_name), to_fn)
  }
  
  # cellbender pdfs
  if (file.exists(sprintf("%s/runs/%s/%s/outs/cellbender_output.pdf", base_dir, run_name, sample_name))) {
    dir.create(sprintf("%s/runs/%s/report/CellBender", base_dir, run_name), showWarnings = F)
    to_fn = sprintf("%s/runs/%s/report/CellBender/%s_cellbender_output.pdf", base_dir, run_name, sample_name)
    if (!file.exists(to_fn)) {
      file.copy(sprintf("%s/runs/%s/%s/outs/cellbender_output.pdf", base_dir, run_name, sample_name), to_fn)
    }
  }
  
  qc_fn = sprintf("%s/runs/%s/report/tmp/%s_qcs.csv", base_dir, run_name, sample_name)
  stopifnot(file.exists(qc_fn))
  qc_m = read.csv(qc_fn, header = T)

  cr_fn = sprintf("%s/runs/%s/%s/outs/metrics_summary.csv", base_dir, run_name, sample_name)
  stopifnot(file.exists(cr_fn))
  cr_m = read.csv(cr_fn, header = T, check.names = F)
  stopifnot(all(selected_cr_columns %in% colnames(cr_m)))
  cr_m = cr_m[, selected_cr_columns]
  perc = grepl("%", cr_m[1,])
  cr_m[, !perc] = as.numeric(gsub(",", "", cr_m[, !perc]))
  cr_m[, perc]  = as.numeric(gsub("%", "", cr_m[, perc])) / 100

  qc_metrics = rbind(qc_metrics, cbind(qc_m, cr_m))
}  

exp_info_fn = sprintf("%s/runs/%s/report/exp_metrics.csv", base_dir, run_name)
stopifnot(file.exists(exp_info_fn))
exp_info = read.csv(exp_info_fn, header=T)
stopifnot(nrow(exp_info) == nrow(qc_metrics) && all(exp_info$Sample.name %in% qc_metrics$Sample.name))

exp_f_cells = c(SC3Pv3 = 0.625, SC3Pv2 = 0.575)
exp_f_doublets = c(SC3Pv3 = 4.8e-6, SC3Pv2 = 4.4e-6)

exp_info = exp_info %>% 
  mutate(Exp_total_cells  = Loaded_total_cells  * exp_f_cells[Chemistry],
         Exp_viable_cells = Loaded_viable_cells * exp_f_cells[Chemistry],
         Exp_f_total_multiplets = Loaded_total_cells * exp_f_doublets[Chemistry],
         Exp_f_viable_multiplets = Loaded_viable_cells * exp_f_doublets[Chemistry]) %>%
  left_join(qc_metrics)

write.csv(exp_info, sprintf("%s/runs/%s/report/full_exp_metrics.csv", base_dir, run_name), row.names=F)

meta_fn = sprintf("%s/scRNAseq_samples_metrics.csv", base_dir)

#lock_fn = sprintf("%s/.scRNAseq_samples_metrics.lock", base_dir)
#st = Sys.time()
#while (file.exists(lock_fn) && Sys.time() - st < 10) {
#  message("Someone else is updating master table, waiting ...")
#  Sys.sleep(2)
#}

#if (file.exists(lock_fn)) {
#  message(sprintf("Error: couldn't lock master table for edit, to force edit delete lock file: %s", lock_fn))
#} else {
#  file.create(lock_fn)
  if (file.exists(meta_fn)) {
    message("Reading master table")
    prev_metrics = read.csv(meta_fn, header=T, check.names = F)
    matched = prev_metrics %>% filter(Run_name == run_name & Sample.name %in% exp_info$Sample.name)
    if (nrow(matched) > 0) {
      message(sprintf("Warning: entries for %d of %d samples already exist in %s/scRNAseq_samples_metrics.csv, not adding any metrics for run %s.\n", nrow(matched), nrow(exp_info), base_dir, run_name))
    } else {
      message("no matched previous samples")
      stopifnot(all(sort(colnames(prev_metrics)) == sort(colnames(exp_info))))
      prev_metrics = rbind(prev_metrics, exp_info[, colnames(prev_metrics)])
      write.csv(prev_metrics, meta_fn, quote=F, row.names=F)
    }
  } else {
    write.csv(exp_info, meta_fn, quote=F, row.names = F)
  }
#  file.remove(lock_fn)
#}
