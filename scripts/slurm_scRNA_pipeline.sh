#!/usr/bin/env bash

#SBATCH -J cc.scRNA
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --signal=2
#SBATCH --mem=128G
#SBATCH -o SCRNASEQ.%a.out
#SBATCH -e SCRNASEQ.%a.err

set -xeu

# Load specific run's parameters
WD=$(pwd -P)
PARAMS_FN=$WD/scRNA_params.sh
if [ ! -e $PARAMS_FN ]
then
        echo "Error: missing parameters file ($PARAMS_FN)"
        exit 1
else
        source $PARAMS_FN
fi

if [ "$WD" != "$BASE_DIR/runs/$RUN_NAME" ]
then
        echo "Error: you should run this from the specific run directory ($BASE_DIR/runs/$RUN_NAME)"
        exit 1
fi

TMP_DIR=/tmp/tmp.${USER}.${SLX_RUN}

META_FILE=$(ls -1 $BASE_DIR/runs/$RUN_NAME/metadata/*csv)
SAMPLE_NAME=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $4 }' | tr -d '"')
SAMPLE_ID=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $2 }' | tr -d '"')

# Run cellranger
if [ "$INCLUDE_INTRONS" == "by_sample_name" ]
then
        if [[ $SAMPLE_NAME =~ "nuclei" ]]
        then
                INCLUDE_INTRONS="--include-introns"
        else
                INCLUDE_INTRONS=""
        fi
fi

if [ "${RUN_CELLRANGER}" == "rerun"  -a -e ${WD}/${SAMPLE_NAME} ]
then
        echo "*** Deleting ${SAMPLE_NAME} folder and reruning cellranger count ***"
        rm -rf ${WD}/${SAMPLE_NAME}
fi

if [ "${RUN_CELLRANGER}" == "rerun" ] || [ "${RUN_CELLRANGER}" == "run" -a ! -e ${WD}/${SAMPLE_NAME}/outs/metrics_summary.csv ]
then
        $BASE_DIR/software/cellranger-5.0.0/cellranger count \
                --id ${SAMPLE_NAME} \
                --fastqs ${WD}/fastq/raw \
                --sample ${SAMPLE_ID} \
                --transcriptome ${REF_CELLRANGER_PATH} \
	        --chemistry $CHEM \
	        $INCLUDE_INTRONS \
                --localcores $SLURM_CPUS_PER_TASK \
                --localmem 128 \
                --maxjobs 80 \
                --jobmode slurm
fi

if [ "${RUN_CELLBENDER}" == "rerun" ] || [ "${RUN_CELLBENDER}" == "run" -a ! -e ${WD}/${SAMPLE_NAME}/outs/cellbender_output_filtered.h5 ]
then
        echo "*** Removing ambient noise with CellBender (sending job to a GPU node) ***"
        ${BASE_DIR}/scripts/scRNA_run_CellBender.R $WD $SAMPLE_NAME \
        $EXPECTED_CELLS $TOTAL_DROPLETS_INCLUDES \
        $FPR $EPOCHS $LEARNING_RATE $SLURM_ARRAY_TASK_ID
fi

# QC filter cells by number of UMIs and %mito UMIs
if [ "${RUN_QC_FILTER}" == "rerun" ] || [ "${RUN_QC_FILTER}" == "run" -a ! -e ${WD}/${SAMPLE_NAME}/outs/qc_filt_barcodes.tsv.gz ]
then
        ${BASE_DIR}/scripts/scRNA_qc_filter_cells.R $WD $SAMPLE_NAME \
        $MIN_HUMAN_NON_MT_UMIS $MAX_HUMAN_F_MITO \
        $MIN_MOUSE_NON_MT_UMIS $MAX_MOUSE_F_MITO \
        $HUMAN_PREF $MOUSE_PREF $MITO_FILTER_BY
fi

# Run souporcell
if [ "${RUN_SOC}" == "rerun"  -a -e ${WD}/${SAMPLE_NAME}/outs/souporcell ]
then
        echo "*** Deleting ${SAMPLE_NAME}/outs/souporcell folder and reruning Souporcell ***"
        rm -rf ${WD}/${SAMPLE_NAME}/outs/souporcell
fi

N_BARCODES=$(zcat ${WD}/${SAMPLE_NAME}/outs/qc_filt_barcodes.tsv.gz | wc -l)
if [ "${RUN_SOC}" == "rerun" ] || [ "${RUN_SOC}" == "run" -a "${SOC_N_CLUSTERS}" -gt 0 -a ! -e ${WD}/${SAMPLE_NAME}/outs/souporcell/consensus.done ]
then
	if (( $N_BARCODES <= 50 ))
	then
		echo "Warning: only ${N_BARCODES} passed QC, too few for souporcell, skipping."
	else	
        	cd $BASE_DIR
		singularity exec -B $BASE_DIR:/dummy_dir $BASE_DIR/software/souporcell.sif souporcell_pipeline.py \
		   -i /dummy_dir/runs/${RUN_NAME}/${SAMPLE_NAME}/outs/possorted_genome_bam.bam \
		   -b /dummy_dir/runs/${RUN_NAME}/${SAMPLE_NAME}/outs/qc_filt_barcodes.tsv.gz \
		   -f /dummy_dir/${SOC_REF_FASTA} \
		   -t $SLURM_CPUS_PER_TASK -o /dummy_dir/runs/${RUN_NAME}/${SAMPLE_NAME}/outs/souporcell -k $SOC_N_CLUSTERS --ignore True
		cd ${WD}
	fi
fi

# Integrate souporcell clustering into the mat and generate some plots
if [ "${RUN_MAT_INTEGRATION}" == "rerun" ] || [ "${RUN_MAT_INTEGRATION}" == "run" -a "${SOC_N_CLUSTERS}" -gt 0 -a ! -e ${WD}/${SAMPLE_NAME}/outs/scrnaseq.done ]
then
        ${BASE_DIR}/scripts/scRNA_SoC_finalise.R $WD $SAMPLE_NAME $HUMAN_PREF $MOUSE_PREF
fi

# Run souporcell on human QC filtered singlet cells (optional)
if [ "${RUN_SOC_ON_HG_SINGLETS}" == "rerun"  -a -e ${WD}/${SAMPLE_NAME}/outs/souporcell_hg ]
then
        echo "*** Deleting ${SAMPLE_NAME}/outs/souporcell_hg folder and reruning Souporcell on human QC filtered singlets ***"
        rm -rf ${WD}/${SAMPLE_NAME}/outs/souporcell_hg
fi

N_HG_BARCODES=$(zcat ${WD}/${SAMPLE_NAME}/outs/qc_hg_singlet_filt_barcodes.tsv.gz | wc -l)
if [ "${RUN_SOC_ON_HG_SINGLETS}" == "rerun" ] || [ "${RUN_SOC_ON_HG_SINGLETS}" == "run" -a ! -e ${WD}/${SAMPLE_NAME}/outs/souporcell_hg/consensus.done ]
then
	if (( $N_HG_BARCODES <= 50 ))
	then
		echo "Warning: only ${N_HG_BARCODES} human singlet cells passed QC, too few for souporcell, skipping."
	else	
    cd $BASE_DIR
		singularity exec -B $BASE_DIR:/dummy_dir $BASE_DIR/software/souporcell.sif souporcell_pipeline.py \
		   -i /dummy_dir/runs/${RUN_NAME}/${SAMPLE_NAME}/outs/possorted_genome_bam.bam \
		   -b /dummy_dir/runs/${RUN_NAME}/${SAMPLE_NAME}/outs/qc_hg_singlet_filt_barcodes.tsv.gz \
		   -f /dummy_dir/${SOC_REF_FASTA} \
		   -t $SLURM_CPUS_PER_TASK -o /dummy_dir/runs/${RUN_NAME}/${SAMPLE_NAME}/outs/souporcell_hg -k $SOC_N_CLUSTERS --ignore True
		cd ${WD}
	fi
fi

