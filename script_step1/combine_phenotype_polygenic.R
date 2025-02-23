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
  FID_temp[i] = str_sub(rownames(polygenic)[i], 1, str_locate(rownames(polygenic)[i], "_")-1)[1]
  IID_temp[i] = str_sub(rownames(polygenic)[i], 1, str_locate(rownames(polygenic)[i], "_")-1)[2]
}


if (all.equal(as.integer(FID_temp), sample$FID) & all.equal(as.integer(IID_temp), sample$IID)) {

    print(head(polygenic))
    print(head(sample))
    print(head(covariate))

} else {

    polygenic = polygenic[match(sample$IID, IID_temp), ]

    FID_temp = IID_temp = c()
    for (i in 1:nrow(polygenic)) {
        FID_temp[i] = str_sub(rownames(polygenic)[i], 1, str_locate(rownames(polygenic)[i], "_")-1)[1]
        IID_temp[i] = str_sub(rownames(polygenic)[i], 1, str_locate(rownames(polygenic)[i], "_")-1)[2]
    }

    if (all.equal(as.integer(FID_temp), sample$FID) & all.equal(as.integer(IID_temp), sample$IID)) {
        print(head(polygenic))
        print(head(sample))
        print(head(covariate))
    }

}

### make sure three data have the same ID

phenotype_all = cbind(sample, polygenic, covariate[, 3:ncol(covariate)])
colnames(phenotype_all)[5:27] = paste0("Chr", 1:23)
head(phenotype_all)

fwrite(phenotype_all, file = "/phenotype_all.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
