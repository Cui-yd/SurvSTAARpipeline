##########################################################
# Split WGS/WES data analysis into specific parts
# Yidan Cui
# Initiate date: 2025/02/21
# Current date: 2025/02/21
##########################################################

rm(list = ls())
gc()


library(SurvSTAAR)
library(data.table)

## if all chromosome are in one file

plink.file = "plink.file.name"
bim_data = fread(paste0(plink.file, ".bim"))

job_num = c()
for (chr in 1:22) {
  chr_num = length(which(bim_data[,1] == chr))
  job_num[chr] = ceiling(chr_num / 1e5)
}

## if chromosome are individually

plink.file = "plink.file.name"  ## without chromosome number

for (chr in 1:22) {
  bim_data_chr = fread(paste0("plink.file"), chr, ".bim")
  chr_num = nrow(bim_data_chr)
  job_num[chr] = ceiling(chr_num / 1e5)
}

job_num
## For example (randomly sampled)
## 24 17 16 19 18 20 14  9 23 14 15 12 12 15 14 12  8  8 13 24 16 17
sum(job_num)
## 340
