#!/usr/bin/env bash

#SBATCH -J cc.scRNA.sw
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --signal=2
#SBATCH --mem=96G
#SBATCH -o SWAPPED_BCS.out
#SBATCH -e SWAPPED_BCS.err

set -xeu

# Load specific run's parameters
WD=$PWD
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

if [ "${RUN_CLEAN_BARCODE_SWAPPING}" == "rerun" ] || [ "${RUN_CLEAN_BARCODE_SWAPPING}" == "run" -a ! -e ${WD}/cleaned_swapped_barcodes.done ]
then
        ${BASE_DIR}/runs/${RUN_NAME}/scRNAseq_pipeline/scripts/scRNA_clean_swapped_barcodes.R ${BASE_DIR}/runs/${RUN_NAME} $SWAPPED_DROPS_MIN_FRAC
fi
