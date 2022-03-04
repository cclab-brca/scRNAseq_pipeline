# Run number (from genomics, must update!)
SLX_RUN=__SLX_RUN__

# Run name (will be the directory name)
RUN_NAME=__RUN_NAME__

BASE_DIR=/mnt/scratchb/cclab/scRNAseq

# Which steps to run (only if their output does not exist), skip, and rerun (run even if output exists. Rerunning cellranger will delete the output folder and will will require rerunning the other steps too. Reruning souporcell will delete the outs/souporcell folder and rerun it from scratch)
RUN_CELLRANGER=run
RUN_CELLBENDER=run
RUN_QC_FILTER=run
RUN_SOC=run
RUN_MAT_INTEGRATION=run

#### Cell ranger reference genome, one of the following:
# refdata-gex-GRCh38-2020-A (human)
# refdata-gex-GRCh38-and-mm10-2020-A (human and mouse)
# refdata-gex-GRCh38-and-mm10-2020-A-eGFP (human, mouse and eGFP)
# refdata-gex-mm10-2020-A (mouse)
REF_CELLRANGER_PATH=$BASE_DIR/ref_data/refdata-gex-GRCh38-and-mm10-2020-A-eGFP

#### 10x library chemistry. Options are: 'auto' for autodetection, 'threeprime' for Single Cell 3', 
# 'fiveprime' for  Single Cell 5', 'SC3Pv1' or 'SC3Pv2' or 'SC3Pv3' for Single Cell 3' v1/v2/v3, 
# 'SC5P-PE' or 'SC5P-R2' for Single Cell 5', paired-end/R2-only, 'SC-FB' for Single Cell Antibody-only 3' v2 or 5' [default: auto] 
CHEM=SC3Pv3

# Count reads mapped to introns? (leave empty for scRNAseq, "--include-introns" if single-nucleus RNA-seq, and "by_sample_name" to use introns if SAMPLE_NAME contains "nuclei")
INCLUDE_INTRONS=""

#### CellBender: remove ambient noise and recall non-empty cells
# Requires installing an anaconda environment named CellBenderGPU (see instructions in the README file).
# By defailt the 'auto' values will use the number of valid cells from cellranger count as the expected cells (and 1.5 that as the total droplets included. make sure this includes empty droplets).
# If you'll set a number in either it will use that number for all samples in the run. 
# For sample specific parameters fill in the parameters in a csv file named reports/CellBender_params.csv with these columns: Sample.name, Expected.cells, Total.droplets.included. You can only fill-in params for some samples, params for the samples not in the csv will be considered as 'auto'. You can also set one of the params and set the other as 'auto'.
EXPECTED_CELLS=auto
TOTAL_DROPLETS_INCLUDES=auto

# Default Cellbender parameters that you'll probably won't need to change
FPR=0.01
EPOCHS=150
LEARNING_RATE=0.0001

#### QC filter cutoffs - by default they are all 'auto' for automatic detection of cutoffs but can be changed to actual values (e.g. MIN_HUMAN_NON_MT_UMIS=750, MAX_HUMAN_F_MITO=0.5)
# To set sample specific values first run the RUN_QC_FILTER step with 'auto', it will create csv files per sample under the report/tmp folder with the automatic detected values.
# Edit these csv files (column names hg/mm_min_nonMT_umi_cutoff and hg/mm_max_f_MT_cutoff), set the below parameters to 'from_table' and rerun the step (set the above RUN_QC_FILTER=rerun)
MIN_HUMAN_NON_MT_UMIS=auto
MIN_MOUSE_NON_MT_UMIS=auto
MAX_HUMAN_F_MITO=auto
MAX_MOUSE_F_MITO=auto

# Filter by mito policy: frac (default, only by max %mito umis), frac_and_count (by max %mito umis and a max count of mito umis, aims to capture doublets of mm and mito-only-hg that are missed by souporcell).
MITO_FILTER_BY=frac

# Gene prefixes (type NA if specie is not expected, e.g. in mouse or human only samples)
HUMAN_PREF=GRCh38
MOUSE_PREF=mm10__

##### Souporcell (SoC)
# SoC reference fasta (expecting local path). Ideally same as REF_CELLRANGER_PATH, but currently SoC cannot handle the combind hg-mm genome, so using only the human genome.
SOC_REF_FASTA=ref_data/refdata-gex-GRCh38-2020-A/fasta/genome.fa
# SoC number of clusters (to skip running SoC enter SOC_N_CLUSTERS=0)
SOC_N_CLUSTERS=2

