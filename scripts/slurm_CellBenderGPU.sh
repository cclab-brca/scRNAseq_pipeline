#!/usr/bin/env bash
#SBATCH --gres gpu:1
#SBATCH -p gpunodes
#SBATCH -J cc.scRNA
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --signal=2
#SBATCH --mem=64G
#SBATCH -o CELLBENDER.%a.out
#SBATCH -e CELLBENDER.%a.err

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

if [ "$WD" != "$BASE_DIR/runs/$RUN_NAME"]
then
        echo "Error: you should run this from the specific run directory ($BASE_DIR/runs/$RUN_NAME)"
        exit 1
fi

META_FILE=$(ls -1 $BASE_DIR/runs/$RUN_NAME/metadata/*csv)
SAMPLE_NAME=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $4 }' | tr -d '"')

CONDA_BIN_DIR=$(dirname $CONDA_EXE)
if [[ -z "$CONDA_BIN_DIR" ]]
then
        echo "Error: Environment variable CONDA_EXE is not defined, make sure you have installed anaconda"
        exit 1
fi

CELLBENDER_ENV="$(dirname $(dirname $CONDA_EXE))/envs/CellBenderGPU"
if [[ ! -d "$CELLBENDER_ENV" ]]
then
        echo "Error: conda environment named CellBenderGPU not found among your conda environments"
        exit 1
fi

set +eu
source ${CONDA_BIN_DIR}/activate CellBenderGPU
set -xeu

cd ${WD}/${SAMPLE_NAME}/outs

cellbender remove-background \
	--input raw_feature_bc_matrix.h5 \
	--output cellbender_output.h5 \
	--cuda \
	--expected-cells $N_CELLS \
	--total-droplets-included $N_DROPLETS \
	--fpr $FPR \
	--epochs $EPOCHS \
	--learning-rate $LEARNING_RATE
