#!/home/bioinformatics/software/R/R-4.0.4-gcc-9.2.0/bin/Rscript
require(dplyr, quietly=T)

# Create a template csv file placeholders for experimental information on the run's samples

args = commandArgs(trailingOnly = T)

base_dir = args[1]
run_name = args[2]

samp_md_fn = list.files(sprintf("%s/runs/%s/metadata", base_dir, run_name), pattern = "*csv", full.names = T)
stopifnot(length(samp_md_fn) == 1)

samp_md = read.csv(samp_md_fn) %>%
  mutate(User = Sys.getenv("USER"),
         Run_name = run_name,
         Treatment = "None",
         Tissue_state = "Frozen",
         Cells_or_Nuclei = "Cells",
         Dissociation = "Miltenyi",
         FACS_sort = "Y",
         Processing_time_in_minutes = "",
         Submission_date = format(Sys.time(), "%d/%m/%Y"),
         cDNA_prep_date = format(Sys.time(), "%d/%m/%Y"),
         Cell_counter = "Sorting",
         Chemistry = "SC3Pv3",
         Loaded_total_cells = "",
         Loaded_viable_cells = "",
         Viability = "",
         Seq_machine = "NovaSeq",
         Flowcell = "",
         Num_lanes = "")

write.csv(samp_md, sprintf("%s/runs/%s/report/exp_metrics.csv", base_dir, run_name), quote = F, row.names = F)
         
