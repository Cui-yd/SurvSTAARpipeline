#!/bin/bash

regenie \
  --step 1 \
  --bed mergeGenotype_prune \
  --covarFile ./path/to/covariates.txt \
  --phenoFile ./path/to/binary_phenotype.txt \
  --bsize 2000 \
  --bt --lowmem \
  --lowmem-prefix tmp_rg \
  --out ./regenie_result