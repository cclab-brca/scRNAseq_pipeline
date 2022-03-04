################################
# Makefile for scRNAseq pipeline
################################
BASE_DIR = /mnt/scratchb/cclab/scRNAseq
TMP_DIR = /tmp/tmp.$(USER).$(SLX_RUN)
CLARITY_JAR = $(BASE_DIR)/software/clarity-tools.jar

# 1. Init
# 
# To run: make init RUN_NAME=<...> SLX_RUN=<...>
# Creates the output folder under runs/$(RUN_NAME) and downloads the metadata file under metadata folder. Check that
# there is only one metadata file there and that it has the correct information.
init:
	mkdir -p runs/$(RUN_NAME)/metadata
	mkdir -p runs/$(RUN_NAME)/fastq/raw
	java -jar $(CLARITY_JAR) -l $(SLX_RUN) -f '*.csv' -d $(TMP_DIR)
	mv $(TMP_DIR)/$(SLX_RUN)/*.csv runs/$(RUN_NAME)/metadata
	ln -sf $(BASE_DIR)/scripts/slurm_scRNA_pipeline.sh runs/$(RUN_NAME)/slurm_scRNA_pipeline.sh
	ln -sf $(BASE_DIR)/scripts/slurm_scRNA_get_fastqs.sh runs/$(RUN_NAME)/slurm_scRNA_get_fastqs.sh
	ln -sf $(BASE_DIR)/scripts/slurm_CellBenderGPU.sh runs/$(RUN_NAME)/slurm_CellBenderGPU.sh
	cat scripts/scRNA_params.sh \
	| sed 's/__RUN_NAME__/$(RUN_NAME)/g' \
	| sed 's/__SLX_RUN__/$(SLX_RUN)/g' \
	> runs/$(RUN_NAME)/scRNA_params.sh
	chmod 755 runs/$(RUN_NAME)/scRNA_params.sh
	make init_report BASE_DIR=$(BASE_DIR) RUN_NAME=$(RUN_NAME)

init_report:
	mkdir -p runs/$(RUN_NAME)/report/tmp
	$(BASE_DIR)/scripts/scRNA_create_run_csv.R $(BASE_DIR) $(RUN_NAME)
	java -jar $(CLARITY_JAR) -l $(SLX_RUN) -f '$(SLX_RUN).*.html' -d $(TMP_DIR)
	mv $(TMP_DIR)/$(SLX_RUN)/$(SLX_RUN).*.html runs/$(RUN_NAME)/report
	
# Throw all qc figs to ppt slides
qcs_to_ppt:
	$(BASE_DIR)/scripts/scRNA_qcs_to_ppt.R $(BASE_DIR)/runs/$(RUN_NAME) $(HG_MM_MODE)
	
# Collect html and csv reports on the run level and all runs meta-table ($(BASE_DIR)/scRNAseq_samples_metrics.csv).
# Entries of the current run will be added to the meta-table if they are not already there (identified by $RUN_NAME and sample names).
report: qcs_to_ppt
	$(BASE_DIR)/scripts/scRNA_make_report.R $(BASE_DIR) $(RUN_NAME)
	
clean_fastq:
	rm -f runs/$(RUN_NAME)/fastq/raw/*fastq.gz
	
clean_all: clean_fastq
	rm -f runs/$(RUN_NAME)/*/outs/*bam runs/$(RUN_NAME)/*/outs/souporcell/souporcell_minimap_tagged_sorted.bam



  

