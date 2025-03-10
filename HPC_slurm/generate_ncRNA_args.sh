## chr_Nsets is the print of "split_jobs_ncRNA.R"
## For example, if you run each job with a maximum of 50 genes:
chr_Nsets=(0 23 20 14 13 16 13 12 13 10 11 14 16  8 11 12 14 16  8 12  8  5  7)

for chr in {1..22}; do
    Nsets=${chr_Nsets[$chr]}
    for set in $(seq 1 $Nsets); do
        echo "${chr} ${set} ${Nsets} ncRNA_results_chr${chr}_${set} objNull_chr${chr}.rda ukb.500k.wgs.chr${chr}.pass.annotated.gds Annotation_name_catalog.rda annotation/info/QC_label TRUE"
    done
done  >> args_ncRNA.txt

