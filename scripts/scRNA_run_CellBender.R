#!/home/bioinformatics/software/R/R-4.0.4-gcc-9.2.0/bin/Rscript

# Wrapper for sending a slurm job on a GPU node for CellBender
require(dplyr, quietly = T)
require(tibble, quietly = T)

args = commandArgs(trailingOnly = T)
if (length(args) == 0) {
  message("Usage: scRNA_run_CellBender.R <base dir> <sample name> <expected cells> <total droplets included> <fpr> <epochs> <learning rate> <sample number>\n")
  stop()
}

base_dir = args[1]
sample_name = args[2]
fpr = as.numeric(args[5])
epochs = as.numeric(args[6])
learning_rate = as.numeric(args[7])
sample_id = args[8]

p_fn = sprintf("%s/report/CellBender_params.csv", base_dir)
if ("from_table" %in% args[3:4]) {
  if (file.exists(p_fn)) {
    ps = read.csv(p_fn) %>% tibble::column_to_rownames('Sample.name')
  } else {
    stop(sprintf("Error: some parameters are required from a sample specific table (%s) which does not exist.", p_fn))
  }
}

samp_dir = sprintf("%s/%s/outs", base_dir, sample_name)

cellranger_fn = sprintf("%s/metrics_summary.csv", samp_dir)
if (!file.exists(cellranger_fn)) {
  stop(sprintf("Error: %s does not exist, make sure you ran cellranger count", cellranger_fn))
}

cr_metrics = read.csv(cellranger_fn, check.names = F)

n_cells = as.numeric(gsub(",", "", cr_metrics[, 'Estimated Number of Cells']))
if (args[3] != 'auto') {
  if (args[3] == 'from_table') {
    if (sample_name %in% rownames(ps) && ps[sample_name, 'Expected.cells'] != 'auto') {
      n_cells = ps[sample_name, 'Expected.cells']
    }
  } else {
    n_cells = as.numeric(args[3])
  }
}

n_droplets = ceiling(n_cells * 1.5)
if (args[4] != 'auto') {
  if (args[4] == 'from_table') {
    if (sample_name %in% rownames(ps) && ps[sample_name, 'Total.droplets.included'] != 'auto') {
      n_droplets = ps[sample_name, 'Total.droplets.included']
    }
  } else {
    n_droplets = as.numeric(args[4])
  }
}

rv = system(sprintf("sbatch -t 6:0:0 -a %s --wait --export=ALL,N_CELLS='%d',N_DROPLETS='%d',FPR='%f',EPOCHS='%d',LEARNING_RATE='%f' slurm_CellBenderGPU.sh", sample_id, n_cells, n_droplets, fpr, epochs, learning_rate))
if (rv > 0) {
  stop("Error: slurm_CellBenderGPU.sh job failed")
}

write.csv(data.frame(Sample.name=sample_name,CellBender_Ncells=n_cells,CellBender_Ndroplets=n_droplets,CellBender_FPR=fpr,CellBender_Epochs=epochs,CellBender_LearnRate=learning_rate), 
          sprintf("%s/report/tmp/%s_qcs.csv", base_dir, sample_name), quote=F, row.names = F)
