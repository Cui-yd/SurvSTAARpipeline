#!/bin/bash

#SBATCH -J indiGDS
#SBATCH -p 64c512g
#SBATCH --ntasks-per-node=20
#SBATCH -o indiGDS_%A_%a.out
#SBATCH -e indiGDS_%A_%a.err


module load miniconda3
source activate survstaar


mapfile -t args < args_IndividualGDS.txt
IFS=' ' read -r chr set Nsets output objnull agds QC_label verbose <<< "${args[$SLURM_ARRAY_TASK_ID-1]}"


Rscript IndividualAnalysisGDS_onCluster.R --chr="$chr" --set="$set" --Nsets="$Nsets" --output.file="$output" --objnull.file="$objnull" --agds.file="$agds" --QC_label="$QC_label" --verbose="$verbose"

