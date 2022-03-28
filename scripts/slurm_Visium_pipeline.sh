#!/usr/bin/env bash

#SBATCH -J visium
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --signal=2
#SBATCH --mem=64G
#SBATCH -o VISIUM.%a.out
#SBATCH -e VISIUM.%a.err

set -xeu

# Load specific run's parameters
WD=$(pwd -P)
PARAMS_FN=$WD/Visium_params.sh
if [ ! -e $PARAMS_FN ]
then
        echo "Error: missing parameters file ($PARAMS_FN)"
        exit 1
else
        source $PARAMS_FN
fi

if [ "$WD" != "$BASE_DIR/runs/Visium/$RUN_NAME" ]
then
  echo "Error: you should run this from the specific run directory ($BASE_DIR/runs/Visium/$RUN_NAME)"
  exit 1
fi

SAMPLE_NAME=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $4 }' | tr -d '"')
SAMPLE_ID=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $2 }' | tr -d '"')

N_MATCHED_SAMPLES_IN_SLIDE_DATA=$(grep $SAMPLE_NAME $SLIDE_DATA_FILE | wc -l)
if [ $N_MATCHED_SAMPLES_IN_SLIDE_DATA -ne 1 ]
then
  echo "Error: multiple lines in slide data file $SLIDE_DATA_FILE match sample $SAMPLE_NAME"
  exit 1
fi

# Individual_image_name field
IMAGE_FILE=$(grep $SAMPLE_NAME $SLIDE_DATA_FILE | awk -F',' '{ print $2 }' | tr -d '"')

# Serial_number field
SN=$(grep $SAMPLE_NAME $SLIDE_DATA_FILE | awk -F',' '{ print $4 }' | tr -d '"')

# Area field
AREA=$(grep $SAMPLE_NAME $SLIDE_DATA_FILE | awk -F',' '{ print $5 }' | tr -d '"')

# Run spaceranger
if [ "${RUN_SPACERANGER}" == "rerun"  -a -e ${WD}/${SAMPLE_NAME} ]
then
        echo "*** Deleting ${SAMPLE_NAME} folder and reruning cellranger count ***"
        rm -rf ${WD}/${SAMPLE_NAME}
fi

if [ "${RUN_SPACERANGER}" == "rerun" ] || [ "${RUN_SPACERANGER}" == "run" -a ! -e ${WD}/${SAMPLE_NAME}/outs/metrics_summary.csv ]
then
$BASE_DIR/software/spaceranger-1.3.1/spaceranger count \
  --id=${SAMPLE_NAME} \
  --fastqs=${WD}/fastq/raw \
  --sample=${SAMPLE_ID} \
  --transcriptome=${REF_SPACERANGER_PATH} \
  --image=${WD}/images/$IMAGE_FILE \
  --slide=$SN \
  --area=$AREA \
  --localcores=$SLURM_CPUS_PER_TASK \
  --localmem=$SLURM_MEM_PER_NODE \
  --maxjobs 20 \
  --jobmode slurm
fi
