##########################################################
# Split WGS/WES data analysis into specific parts
# Yidan Cui
# Initiate date: 2025/02/21
# Current date: 2025/02/21
##########################################################

rm(list = ls())
gc()


library(SurvSTAAR)
library(SeqArray)


aGDS.file = "aGDS.file.name"  ## without chromosome number

for (chr in 1:22) {
  genofile = seqOpen(paste0("aGDS.file"), chr, ".gds")
  variant.id = seqGetData(genofile, "variant.id")
  chr_num = length(variant.id)
  job_num[chr] = ceiling(chr_num / 5e5)
}

job_num
## For example (randomly sampled)
## 24 17 16 19 18 20 14  9 23 14 15 12 12 15 14 12  8  8 13 24 16 17
sum(job_num)
## 340
