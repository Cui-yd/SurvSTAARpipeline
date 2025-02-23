#!/bin/bash

./plink2 \
    --bfile mergeGenotype  \
    --maf 0.01 \
    --geno 0.01 \
    --hwe 1e-15 \
    --indep-pairwise 1000 100 0.9 \
    --out mergeGenotype_prune
