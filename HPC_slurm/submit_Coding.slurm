#!/bin/bash

#SBATCH -J coding
#SBATCH -p 64c512g
#SBATCH --ntasks-per-node=20
#SBATCH -o coding_%A_%a.out
#SBATCH -e coding_%A_%a.err


module load miniconda3
source activate survstaar


mapfile -t args < args_Coding.txt
IFS=' ' read -r chr set Nsets output objnull agds categories Annotation_catalog QC_label verbose <<< "${args[$SLURM_ARRAY_TASK_ID-1]}"


Rscript GeneCentricCoding_onCluster.R --chr="$chr" --set="$set" --Nsets="$Nsets" --output.file="$output" --objnull.file="$objnull" --agds.file="$agds" --categories="$categories" --Annotation.name.catalog.file="$Annotation_catalog" --QC_label="$QC_label" --verbose="$verbose"

