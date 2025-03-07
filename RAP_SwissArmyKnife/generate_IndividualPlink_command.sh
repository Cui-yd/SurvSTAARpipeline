## chr_Nsets is the print of "split_jobs_IndividualPlink.R"
## For example, if you have the following print:
chr_Nsets=(0 24 17 16 19 18 20 14  9 23 14 15 12 12 15 14 12  8  8 13 24 16 17)

for chr in {1..22}; do
    Nsets=${chr_Nsets[$chr]}
    for set in $(seq 1 $Nsets); do
        echo "chr${chr}_set${set}"
    done

    for set in $(seq 1 $Nsets); do
        echo "Rscript IndividualAnalysisPlink_onRAP.R --chr=${chr} --set=${set} --Nsets=${Nsets} --output.file=Individual_results_chr${chr}_${set} --objnull.file=objNull_chr${chr}.rda --plink.file=ukb.plink.chr${chr} --verbose=T"
    done
done  >> SAK_batch_individualPlink.txt
