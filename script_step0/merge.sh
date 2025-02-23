#!/bin/bash

./plink1.9/plink \
    --bfile /ukb22418_c1_b0_v2 \
    --merge-list ./merge_list.txt \
    --threads 40 \
    --make-bed \
    --out ./mergeGenotype

## if you want select a subset individuals, add --keep subset_sample.fam \ in the above script.

