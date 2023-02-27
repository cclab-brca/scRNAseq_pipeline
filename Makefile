################################
# Makefile for scRNAseq pipeline
################################
#BASE_DIR = /mnt/scratcha/jblab/eyal-l01/scRNAseq
BASE_DIR = /mnt/scratchb/cclab/scRNAseq
TMP_DIR = /tmp/tmp.$(USER).$(SLX_RUN)
JAVA=/home/bioinformatics/software/java/java17/bin/java
CLARITY_JAR = $(BASE_DIR)/software/clarity-tools.jar

#############################
# 
# scRNA-seq: using cellranger
# 
#############################
# 1. Init
# 
# To run: make init RUN_NAME=<...> SLX_RUN=<...>
# Creates the output folder under runs/$(RUN_NAME) and downloads the metadata file under metadata folder. Check that
# there is only one metadata file there and that it has the correct information.
init:
	mkdir -p runs/$(RUN_NAME)/metadata
	mkdir -p runs/$(RUN_NAME)/fastq/raw
	$(JAVA) -jar $(CLARITY_JAR) -l $(SLX_RUN) -f '*.csv' -d $(TMP_DIR)
	mv $(TMP_DIR)/$(SLX_RUN)/*.csv runs/$(RUN_NAME)/metadata
	git clone git@github.com:cclab-brca/scRNAseq_pipeline.git runs/$(RUN_NAME)/scRNAseq_pipeline
	for fn in slurm_scRNA_pipeline.sh slurm_scRNA_get_fastqs.sh slurm_CellBenderGPU.sh slurm_scRNA_clean_swapped_barcodes.sh; do \
			ln -s $(BASE_DIR)/runs/$(RUN_NAME)/scRNAseq_pipeline/scripts/$$fn runs/$(RUN_NAME)/$$fn; \
	done
	cat $(BASE_DIR)/scRNAseq_pipeline/scripts/scRNA_params.sh \
	| sed 's/__RUN_NAME__/$(RUN_NAME)/g' \
	| sed 's/__SLX_RUN__/$(SLX_RUN)/g' \
	> runs/$(RUN_NAME)/scRNA_params.sh
	chmod 755 runs/$(RUN_NAME)/scRNA_params.sh
	make init_report BASE_DIR=$(BASE_DIR) RUN_NAME=$(RUN_NAME)

init_report:
	mkdir -p runs/$(RUN_NAME)/report/tmp
	$(BASE_DIR)/runs/$(RUN_NAME)/scRNAseq_pipeline/scripts/scRNA_create_run_csv.R $(BASE_DIR) $(RUN_NAME)
	sed 's/"//g' -i $(BASE_DIR)/runs/$(RUN_NAME)/metadata/*.csv
	$(JAVA) -jar $(CLARITY_JAR) -l $(SLX_RUN) -f '$(SLX_RUN).*.html' -d $(TMP_DIR)
	mv $(TMP_DIR)/$(SLX_RUN)/$(SLX_RUN).*.html runs/$(RUN_NAME)/report
	
# Throw all qc figs to ppt slides
qcs_to_ppt:
	$(BASE_DIR)/runs/$(RUN_NAME)/scRNAseq_pipeline/scripts/scRNA_qcs_to_ppt.R $(BASE_DIR)/runs/$(RUN_NAME) $(HG_MM_MODE)
	
# Collect html and csv reports on the run level and all runs meta-table ($(BASE_DIR)/scRNAseq_samples_metrics.csv).
# Entries of the current run will be added to the meta-table if they are not already there (identified by $RUN_NAME and sample names).
report: qcs_to_ppt
	$(BASE_DIR)/runs/$(RUN_NAME)/scRNAseq_pipeline/scripts/scRNA_make_report.R $(BASE_DIR) $(RUN_NAME)

# Maintenance	
clean_fastq:
	rm -f runs/$(RUN_NAME)/fastq/raw/*fastq.gz
	
clean_all: clean_fastq
	rm -rf runs/$(RUN_NAME)/*/outs/*bam runs/$(RUN_NAME)/*/outs/souporcell/souporcell_minimap_tagged_sorted.bam runs/$(RUN_NAME)/*/SC_RNA_COUNTER_CS

# Copy bams from run dir to a temporary storage (expects RUN_NAME and DEST_DIR to be specified) and vice-versa
copy_bams:
	$(eval SAMPLES := $(shell tail -n +2 runs/$(RUN_NAME)/metadata/*.csv | cut -f 4 -d ','))
	for samp in $(SAMPLES); do \
		echo $$samp; \
		mkdir -p $(DEST_DIR)/$(RUN_NAME)/$$samp/outs/souporcell; \
		cp runs/$(RUN_NAME)/$$samp/outs/possorted_genome_bam.bam $(DEST_DIR)/$(RUN_NAME)/$$samp/outs/ ; \
		cp runs/$(RUN_NAME)/$$samp/outs/souporcell/souporcell_minimap_tagged_sorted.bam $(DEST_DIR)/$(RUN_NAME)/$$samp/outs/souporcell/ ; \
	done
	
copy_back_bams:
	$(eval SAMPLES := $(shell tail -n +2 runs/$(RUN_NAME)/metadata/*.csv | cut -f 4 -d ','))
	for samp in $(SAMPLES); do \
		echo $$samp; \
		cp $(DEST_DIR)/$(RUN_NAME)/$$samp/outs/possorted_genome_bam.bam runs/$(RUN_NAME)/$$samp/outs/  ; \
		cp $(DEST_DIR)/$(RUN_NAME)/$$samp/outs/souporcell/souporcell_minimap_tagged_sorted.bam runs/$(RUN_NAME)/$$samp/outs/souporcell/ ; \
	done
	
	
#############################
# 
# Visium: using spaceranger
# 
#############################
# 1. Init
# 
# To run: make visium_init RUN_NAME=<...> SLX_RUN=<...>
# Creates the output folder under runs/Visium/$(RUN_NAME) and downloads the metadata file under metadata folder. Check that
# there is only one metadata file there and that it has the correct information.
visium_init:
	mkdir -p runs/Visium/$(RUN_NAME)/metadata
	mkdir -p runs/Visium/$(RUN_NAME)/fastq/raw
	$(JAVA) -jar $(CLARITY_JAR) -l $(SLX_RUN) -f '*.csv' -d $(TMP_DIR)
	mv $(TMP_DIR)/$(SLX_RUN)/*.csv runs/Visium/$(RUN_NAME)/metadata
	git clone git@github.com:cclab-brca/scRNAseq_pipeline.git runs/Visium/$(RUN_NAME)/scRNAseq_pipeline
	ln -sf $(BASE_DIR)/runs/$(RUN_NAME)/scRNAseq_pipeline/scripts/slurm_Visium_pipeline.sh runs/Visium/$(RUN_NAME)/slurm_Visium_pipeline.sh
	ln -sf $(BASE_DIR)/runs/$(RUN_NAME)/scRNAseq_pipeline/scripts/slurm_Visium_get_fastqs.sh runs/Visium/$(RUN_NAME)/slurm_Visium_get_fastqs.sh
	cat $(BASE_DIR)/scRNAseq_pipeline/scripts/Visium_params.sh \
	| sed 's/__RUN_NAME__/$(RUN_NAME)/g' \
	| sed 's/__SLX_RUN__/$(SLX_RUN)/g' \
	> runs/Visium/$(RUN_NAME)/Visium_params.sh
	chmod 755 runs/Visium/$(RUN_NAME)/Visium_params.sh
	mkdir -p runs/Visium/$(RUN_NAME)/report/tmp
	$(JAVA) -jar $(CLARITY_JAR) -l $(SLX_RUN) -f '$(SLX_RUN).*.html' -d $(TMP_DIR)
	mv $(TMP_DIR)/$(SLX_RUN)/$(SLX_RUN).*.html runs/Visium/$(RUN_NAME)/report
	

  

