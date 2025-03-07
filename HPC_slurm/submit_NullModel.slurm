#!/bin/bash

#SBATCH -J nullmodel
#SBATCH -p 64c512g
#SBATCH --ntasks-per-node=20
#SBATCH -o nullmodel_%A_%a.out
#SBATCH -e nullmodel_%A_%a.err


module load miniconda3
source activate survstaar


mapfile -t args < args_NullModel.txt
IFS=' ' read -r chr phenofile covCol output verbose <<< "${args[$SLURM_ARRAY_TASK_ID-1]}"


Rscript fitNullModel_onCluster.R --chr="$chr" --phenofile="$phenofile" --covCol="$covCol" --output="$output" --verbose="$verbose"

