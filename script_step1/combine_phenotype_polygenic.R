library(data.table)
library(stringr)

polygenic = fread("./path/to/regenie_result_1.loco")
polygenic = t(polygenic)[-1, ]

sample = fread("./path/to/phenotype.txt")
sample = sample[order(sample$IID),]

covariate = fread("./path/to/covariates.txt")
covariate = covariate[order(covariate$IID),]

FID_temp = IID_temp = c()
for (i in 1:nrow(polygenic)) {
  split_name <- strsplit(rownames(polygenic)[i], "_")[[1]]
  FID_temp[i] <- split_name[1]
  IID_temp[i] <- split_name[2]
}

# Convert to numeric only if you're sure they're all numeric
if (is.numeric(sample$FID) && is.numeric(sample$IID)) {
  FID_temp <- as.numeric(FID_temp)
  IID_temp <- as.numeric(IID_temp)
}

# Check if they're identical (not just equal)
if (identical(FID_temp, sample$FID) && identical(IID_temp, sample$IID)) {
  print(head(polygenic))
  print(head(sample))
  print(head(covariate))
} else {
  # Your reordering code
  polygenic <- polygenic[match(sample$IID, IID_temp), ]

  # Redo FID/IID extraction with corrected logic
  FID_temp = IID_temp = c()
  for (i in 1:nrow(polygenic)) {
    split_name <- strsplit(rownames(polygenic)[i], "_")[[1]]
    FID_temp[i] <- split_name[1]
    IID_temp[i] <- split_name[2]
  }

  # Same check as before
  if (identical(FID_temp, sample$FID) && identical(IID_temp, sample$IID)) {
    print(head(polygenic))
    print(head(sample))
    print(head(covariate))
  }
}

### make sure three data have the same ID

phenotype_all = cbind(sample, covariate[, 3:ncol(covariate)], polygenic)
colnames(phenotype_all)[(ncol(sample)+ncol(covariate)-1) : ncol(phenotype_all)] = paste0("Chr", 1:23)
head(phenotype_all)

fwrite(phenotype_all, file = "/phenotype_all.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
