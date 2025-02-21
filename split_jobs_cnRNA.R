##########################################################
# Split WGS/WES data analysis into specific parts
# Yidan Cui
# Initiate date: 2025/02/20
# Current date: 2025/02/20
##########################################################

rm(list = ls())
gc()


library(SurvSTAAR)

job_num = c()
for (chr in 1:22) {

  ## For ncRNA analysis, please use the following code:
  chr_num = length(which(ncRNA_info$Chr == chr))
  job_num[chr] = ceiling(chr_num / 80)

}

job_num
## For example, if you run each job with a maximum of 80 genes:
## 23 20 14 13 16 13 12 13 10 11 14 16  8 11 12 14 16  8 12  8  5  7
sum(job_num)
## 276
