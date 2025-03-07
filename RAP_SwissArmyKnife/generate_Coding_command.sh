## chr_Nsets is the print of "split_jobs_coding_noncoding.R"
## For example, if you run each job with a maximum of 50 genes:
chr_Nsets=(0 40 25 21 15 18 21 18 14 16 15 26 20  7 12 12 17 23  6 28 11  5  9)

for chr in {1..22}; do
    Nsets=${chr_Nsets[$chr]}
    for set in $(seq 1 $Nsets); do
        echo "chr${chr}_set${set}"
    done

    for set in $(seq 1 $Nsets); do
        echo "Rscript GeneCentricCoding_onRAP.R --chr=${chr} --set=${set} --Nsets=${Nsets} --output.file=Coding_results_chr${chr}_${set} --objnull.file=objNull_chr${chr}.rda --agds.file=ukb.500k.wgs.chr${chr}.pass.annotated.gds --categories=all --Annotation.name.catalog.file=Annotation_name_catalog.rda --QC_label=annotation/info/QC_label --verbose=T"
    done
done  >> SAK_batch_coding.txt
