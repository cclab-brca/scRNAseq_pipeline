# Run number (from genomics, must update!)
SLX_RUN=__SLX_RUN__

# Run name (will be the directory name)
RUN_NAME=__RUN_NAME__

BASE_DIR=/mnt/scratcha/cclab/scRNAseq

# Run metadata files
META_FILE=$BASE_DIR/runs/Visium/$RUN_NAME/metadata/${SLX_RUN}*.contents.csv
SLIDE_DATA_FILE=$BASE_DIR/runs/Visium/$RUN_NAME/metadata/${SLX_RUN}_slide_data.csv

# Which steps to run (only if their output does not exist), skip, and rerun (run even if output exists. Rerunning cellranger will delete the output folder and will will require rerunning the other steps too. Reruning souporcell will delete the outs/souporcell folder and rerun it from scratch)
RUN_SPACERANGER=run

#### spaceranger reference genome, one of the following:
# refdata-gex-GRCh38-2020-A (human)
# refdata-gex-GRCh38-and-mm10-2020-A (human and mouse)
# refdata-gex-GRCh38-and-mm10-2020-A-eGFP (human, mouse and eGFP)
# refdata-gex-mm10-2020-A (mouse)
REF_SPACERANGER_PATH=$BASE_DIR/ref_data/refdata-gex-GRCh38-and-mm10-2020-A-eGFP


#### QC filter cutoffs - by default they are all 'auto' for automatic detection of cutoffs but can be changed to actual values (e.g. MIN_HUMAN_NON_MT_UMIS=750, MAX_HUMAN_F_MITO=0.5)
# To set sample specific values first run the RUN_QC_FILTER step with 'auto', it will create csv files per sample under the report/tmp folder with the automatic detected values.
# Edit these csv files (column names hg/mm_min_nonMT_umi_cutoff and hg/mm_max_f_MT_cutoff), set the below parameters to 'from_table' and rerun the step (set the above RUN_QC_FILTER=rerun)
MIN_HUMAN_NON_MT_UMIS=auto
MIN_MOUSE_NON_MT_UMIS=auto
MAX_HUMAN_F_MITO=auto
MAX_MOUSE_F_MITO=auto


# Gene prefixes (type NA if specie is not expected, e.g. in mouse or human only samples)
HUMAN_PREF=GRCh38
MOUSE_PREF=mm10__


