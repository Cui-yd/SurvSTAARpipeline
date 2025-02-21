##########################################################
# Analyze ncRNA functional categories using SurvSTAAR
# on UK Biobank Research Analysis Platform
# Yidan Cui
# Initiate date: 2025/02/20
# Current date: 2025/02/20
##########################################################

rm(list = ls())
gc()

# set parameters ----------------------------------------------------------

default_args <- list(
  chr = 666,
  set = 666,
  Nsets = 666,
  output.file = "output_default",
  objnull.file = "",
  agds.file = "",
  QC_label = "annotation/filter",
  Annotation.name.catalog.file = "",
  Annotation_dir = "annotation/info/FunctionalAnnotation",
  Use_annotation_weights = TRUE,
  Annotation_name = c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive",
                      "aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                      "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability",
                      "aPC.TF","aPC.Protein"),
  variant_type = "SNV",
  geno_missing_cutoff = 1,
  geno_missing_imputation = "mean",
  rare_maf_cutoff = 0.01,
  rare_num_cutoff = 2,
  min_maf_cutoff = 0,
  combine_ultra_rare = TRUE,
  ultra_rare_mac_cutoff = 20,
  use_SPA = TRUE,
  SPA_filter = TRUE,
  SPA_filter_cutoff = 0.05,
  rm_long = TRUE,
  rm_long_cutoff = 5000,
  verbose = TRUE
)

num_args = c("chr", "set", "Nsets", "geno_missing_cutoff", "rare_maf_cutoff", "rare_num_cutoff",
             "min_maf_cutoff", "ultra_rare_mac_cutoff", "SPA_filter_cutoff", "rm_long_cutoff")
log_args = c("Use_annotation_weights", "combine_ultra_rare", "use_SPA", "SPA_filter", "rm_long", "verbose")


args <- commandArgs(TRUE)

args_list = argsReshape(default_args, args, num_args, log_args)

print(args_list)


print("** Loading packages and data... **")
suppressPackageStartupMessages(library(gdsfmt, quietly = T))
suppressPackageStartupMessages(library(SeqArray, quietly = T))
suppressPackageStartupMessages(library(SeqVarTools, quietly = T))
suppressPackageStartupMessages(library(Matrix, quietly = T))
suppressPackageStartupMessages(library(stats, quietly = T))
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(CompQuadForm, quietly = T))
suppressPackageStartupMessages(library(SurvSTAAR, quietly = T))


print(sessionInfo())


genofile = seqOpen(args_list$agds.file)
print(paste0("Open aGDS file successfully! ", basename(genofile$filename)))

rap_load_as(args_list$objnull.file, "objNull")
print(paste0("Load objnull data successfully! ", args_list$objnull.file))

rap_load_as(args_list$Annotation.name.catalog.file, "Annotation_name_catalog")
print(paste0("Load annotation catalog data successfully! ", args_list$Annotation.name.catalog.file))


ncRNA_info_chr = ncRNA_info[ncRNA_info$Chr == args_list$chr, ]
print(paste0("Number of masks in chr", args_list$chr, " is: ", nrow(ncRNA_info_chr)))


Nsets = args_list$Nsets
set = args_list$set
length_allSet = nrow(ncRNA_info_chr)


if (Nsets != 1) {
  print(paste0("All the masks will be devided into ", args_list$Nsets, " groups"))

  if (set == Nsets) {
    numSet = length_allSet - ceiling(length_allSet / Nsets) * (Nsets - 1)
    idSet = ncRNA_info_chr$ncRNA[(length_allSet - numSet + 1):length_allSet]
  } else {
    numSet = ceiling(length_allSet / Nsets)
    idSet = ncRNA_info_chr$ncRNA[(numSet * (set - 1) + 1):(numSet * set)]
  }

  print(paste0("In sub-job ", set , ", there are ", numSet, " masks will be tested"))

} else if (Nsets == 1) {

  numSet = length_allSet
  print(paste0("In this job, there are ", numSet, " masks will be tested"))

}


print("** All the data and functions required are settled down! **")



# ncRNA analysis --------------------------------------------

output_name = paste0(args_list$output.file,".rda")

ncRNA_results = list()
time_elapse = c()
begin = Sys.time()

print(paste0("Begin running the program at ", begin))

for (j in 1:length(idSet)) {

  begin_j = Sys.time()

  gene_name = idSet[j]
  print(paste0(j, "/", length(idSet), "      ", gene_name))

  ncRNA_results[[j]] = try(ncRNA(gene_name = gene_name,
                                 chr = args_list$chr,
                                 genofile = genofile,
                                 objNull = objNull,
                                 variant_type = args_list$variant_type,
                                 rare_maf_cutoff = args_list$rare_maf_cutoff,
                                 rare_num_cutoff = args_list$rare_num_cutoff,
                                 geno_missing_cutoff = args_list$geno_missing_cutoff,
                                 geno_missing_imputation = args_list$geno_missing_imputation,
                                 min_maf_cutoff = args_list$min_maf_cutoff,
                                 combine_ultra_rare = args_list$combine_ultra_rare,
                                 ultra_rare_mac_cutoff = args_list$ultra_rare_mac_cutoff,
                                 QC_label = args_list$QC_label,
                                 Annotation_dir = args_list$Annotation_dir,
                                 Annotation_name_catalog = Annotation_name_catalog,
                                 Use_annotation_weights = args_list$Use_annotation_weights,
                                 Annotation_name = args_list$Annotation_name,
                                 use_SPA = args_list$use_SPA,
                                 SPA_filter = args_list$SPA_filter,
                                 SPA_filter_cutoff = args_list$SPA_filter_cutoff,
                                 rm_long = args_list$rm_long,
                                 rm_long_cutoff = args_list$rm_long_cutoff,
                                 verbose = args_list$verbose), silent = FALSE)

  time_elapse[j] = difftime(Sys.time(), begin_j, units = "mins")
  print(paste0("Time elapse: ", round(time_elapse[j], 2), " mins"))

  seqResetFilter(genofile)
  gc()

  save(ncRNA_results, file = output_name)

}

end = Sys.time()
print(paste0("Finish analyzing this job at ", end))
elapse = difftime(end, begin, units = "hours")
print(paste0("Time elapse: ", round(elapse, 2), " hours"))
print(paste0("Average time consuming per mask is ", round(mean(time_elapse), 2), " mins"))

seqClose(genofile)
gc()

print(paste0("The results has been saved at ", output_name))


system(paste0("dx upload ", output_name))
