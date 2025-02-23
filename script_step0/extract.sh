#!/bin/bash

./plink2 \
    --bfile mergeGenotype \
    --extract mergeGenotype_prune.prune.in \
    --make-bed \
    --out mergeGenotype_prune
