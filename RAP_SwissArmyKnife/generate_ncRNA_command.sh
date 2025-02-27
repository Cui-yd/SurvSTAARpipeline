## chr_Nsets is the print of "split_jobs_cnRNA.R"
## For example, if you run each job with a maximum of 50 genes:
chr_Nsets=(23 20 14 13 16 13 12 13 10 11 14 16  8 11 12 14 16  8 12  8  5  7)

for chr in {1..22}; do
    Nsets=${chr_Nsets[$chr]}
    for set in $(seq 1 $Nsets); do
        echo "chr${chr}_set${set}"
    done

    for set in $(seq 1 $Nsets); do
        echo "Rscript ncRNA_onRAP.R --chr=${chr} --set=${set} --Nsets=${Nsets} --output.file=ncRNA_results_chr${chr}_${set} --objnull.file=objNull_chr${chr}.rda --agds.file=ukb.500k.wgs.chr${chr}.pass.annotated.gds --Annotation.name.catalog.file=Annotation_name_catalog.rda --QC_label=annotation/info/QC_label --verbose=T"
    done
done  >> SAK_batch_ncRNA.txt
