## chr_Nsets is the print of "split_jobs_coding_noncoding.R"
## For example, if you run each job with a maximum of 50 genes:
chr_Nsets=(0 40 25 21 15 18 21 18 14 16 15 26 20  7 12 12 17 23  6 28 11  5  9)

for chr in {1..22}; do
    Nsets=${chr_Nsets[$chr]}
    for set in $(seq 1 $Nsets); do
        echo "${chr} ${set} ${Nsets} Coding_results_chr${chr}_${set} objNull_chr${chr}.rda ukb.500k.wgs.chr${chr}.pass.annotated.gds all Annotation_name_catalog.rda annotation/info/QC_label TRUE"
    done
done  >> args_Coding.txt

