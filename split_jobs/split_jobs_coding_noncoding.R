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

  ## For gene-centric coding and noncoding analysis
  chr_num = length(which(genes_info$Chr == chr))
  job_num[chr] = ceiling(chr_num / 50)

}

job_num
## For example, if you run each job with a maximum of 50 genes:
## 40 25 21 15 18 21 18 14 16 15 26 20  7 12 12 17 23  6 28 11  5  9
sum(job_num)
## 379
