#!/bin/bash

#SBATCH -J indiPlink
#SBATCH -p 64c512g
#SBATCH --ntasks-per-node=20
#SBATCH -o indiPlink_%A_%a.out
#SBATCH -e indiPlink_%A_%a.err


module load miniconda3
source activate survstaar


mapfile -t args < args_IndividualPlink.txt
IFS=' ' read -r chr set Nsets output objnull plink verbose <<< "${args[$SLURM_ARRAY_TASK_ID-1]}"


Rscript IndividualAnalysisPlink_onCluster.R --chr="$chr" --set="$set" --Nsets="$Nsets" --output.file="$output" --objnull.file="$objnull" --plink.file="$plink" --verbose="$verbose"

