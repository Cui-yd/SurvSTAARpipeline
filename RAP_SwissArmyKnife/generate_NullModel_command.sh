## If you are using the LOCO procedure, please follow this code

for chr in {1..22}; do echo "chr${chr}" >> SAK_batch_null_model.txt; done
for chr in {1..22}; do echo "Rscript fitNullModel_onRAP.R --chr=${chr} --phenofile=phenotype_all.txt  --covCol=sex,birthyr,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10  --output=objNull_chr${chr}  --verbose=T" >> SAK_batch_null_model.txt; done
