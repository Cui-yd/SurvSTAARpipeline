## chr_Nsets is the print of "split_jobs_IndividualPlink.R"
## For example, if you have the following print:
chr_Nsets=(0 24 17 16 19 18 20 14  9 23 14 15 12 12 15 14 12  8  8 13 24 16 17)

for chr in {1..22}; do
    Nsets=${chr_Nsets[$chr]}
    for set in $(seq 1 $Nsets); do
        echo "${chr} ${set} ${Nsets} Individual_results_chr${chr}_${set} objNull_chr${chr}.rda ukb.500k.wgs.chr${chr}.pass.annotated.gds annotation/info/QC_label TRUE"
    done
done  >> args_IndividualGDS.txt

