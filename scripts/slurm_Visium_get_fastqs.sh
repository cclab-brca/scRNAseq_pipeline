#!/usr/bin/env bash

#SBATCH -J visum.fq
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --signal=2
#SBATCH --mem=8G
#SBATCH -o GET_FQS.out
#SBATCH -e GET_FQS.err

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

if [ "$WD" != "$BASE_DIR/runs/Visium/$RUN_NAME" ]
then
        echo "Error: you should run this from the specific run directory ($BASE_DIR/runs/Visium/$RUN_NAME)"
        exit 1
fi

# Initial download to determine input type (old: tar per sample, new: fqs that need renaming)
java -jar ${BASE_DIR}/software/clarity-tools.jar -l $SLX_RUN -f "${SLX_RUN}.*.tar" -f "${SLX_RUN}*.fq.gz" -d ${WD}/fastq/tmp 
if ls ${WD}/fastq/tmp/${SLX_RUN}/*lostreads* >& /dev/null
then
        rm ${WD}/fastq/tmp/${SLX_RUN}/*lostreads*
fi

if ls ${WD}/fastq/tmp/${SLX_RUN}/${SLX_RUN}.*tar >& /dev/null
 then
        mv ${WD}/fastq/tmp/${SLX_RUN}/${SLX_RUN}.*tar $WD/fastq/raw
        for i in $WD/fastq/raw/${SLX_RUN}.*.tar 
        do
                tar xvf $i -C $WD/fastq/raw
                rm $i
        done
else
        mv ${WD}/fastq/tmp/${SLX_RUN}/${SLX_RUN}.*.fq.gz $WD/fastq/raw
        ${BASE_DIR}/scripts/crukci_to_illumina.py ${WD}/fastq/raw
fi

