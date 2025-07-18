library(data.table)
library(stringr)

polygenic = fread("./path/to/regenie_result_1.loco")
polygenic = t(polygenic)[-1, ]

sample = fread("./path/to/phenotype.txt")
sample = sample[order(sample$IID),]

covariate = fread("./path/to/covariates.txt")
covariate = covariate[order(covariate$IID),]



################################################################
## If your ID does not contain underscores "_", use this:
FID_temp = IID_temp = c()
for (i in 1:nrow(polygenic)) {
  split_name <- strsplit(rownames(polygenic)[i], "_")[[1]]
  FID_temp[i] <- split_name[1]
  IID_temp[i] <- split_name[2]
}
## If your ID has a fixed length, use this:
FID_temp = IID_temp = c()
ID_length = nchar(sample$IID[1])
for (i in 1:nrow(polygenic)) {
  split_name <- strsplit(rownames(polygenic)[i], "_")[[1]]
  FID_temp[i] <- substr(rownames(polygenic)[i], 1, ID_length)
  IID_temp[i] <- substr(rownames(polygenic)[i], ID_length + 2, nchar(rownames(polygenic)[i]))
}
## If neither of the above applies, you may need to write custom code to split FID and IID
################################################################



## Check if they're identical
if (identical(FID_temp, sample$FID) && identical(IID_temp, sample$IID)) {
  print(head(polygenic))
  print(head(sample))
  print(head(covariate))
} else {
  polygenic <- polygenic[match(sample$IID, IID_temp), ]
}



###################### check again #############################
## If your ID does not contain underscores "_", use this:
FID_temp = IID_temp = c()
for (i in 1:nrow(polygenic)) {
  split_name <- strsplit(rownames(polygenic)[i], "_")[[1]]
  FID_temp[i] <- split_name[1]
  IID_temp[i] <- split_name[2]
}
## If your ID has a fixed length, use this:
FID_temp = IID_temp = c()
ID_length = nchar(sample$IID[1])
for (i in 1:nrow(polygenic)) {
  split_name <- strsplit(rownames(polygenic)[i], "_")[[1]]
  FID_temp[i] <- substr(rownames(polygenic)[i], 1, ID_length)
  IID_temp[i] <- substr(rownames(polygenic)[i], ID_length + 2, nchar(rownames(polygenic)[i]))
}
## If neither of the above applies, you may need to write custom code to split FID and IID
################################################################



## Final check to ensure the polygenic data and phenotype data are correctly matched
## Make sure three data have the same FID and IID
if (identical(FID_temp, sample$FID) && identical(IID_temp, sample$IID)) {
  print(head(polygenic))
  print(head(sample))
  print(head(covariate))
}


phenotype_all = cbind(sample, covariate[, 3:ncol(covariate)], polygenic)
colnames(phenotype_all)[(ncol(sample)+ncol(covariate)-1) : ncol(phenotype_all)] = paste0("Chr", 1:23)
head(phenotype_all)

fwrite(phenotype_all, file = "/phenotype_all.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
