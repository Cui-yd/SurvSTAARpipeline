##########################################################
# Individual analysis using SurvSTAAR in plink format data
# on High Performance Computing Cluster
# Yidan Cui
# Initiate date: 2025/02/21
# Current date: 2025/02/21
##########################################################

rm(list = ls())
gc()

# set parameters ----------------------------------------------------------

default_args = list(
  chr = 666,
  set = NULL,
  Nsets = NULL,
  output.file = "output_default",
  objnull.file = "",
  plink.file = "",
  start_loc = NULL,
  end_loc = NULL,
  sampleCol = "IID",
  geno_missing_cutoff = 1,
  geno_missing_imputation = "mean",
  use_SPA = TRUE,
  SPA_filter = TRUE,
  SPA_filter_cutoff = 0.05,
  min_mac_cutoff = NULL,
  min_maf_cutoff = 0.01,
  chunk_size = 1e3,
  verbose = TRUE
)

num_args = c("chr", "set", "Nsets", "start_loc", "end_loc", "geno_missing_cutoff",
             "SPA_filter_cutoff", "min_mac_cutoff", "min_maf_cutoff", "chunk_size")
log_args = c("use_SPA", "SPA_filter", "verbose")

args <- commandArgs(TRUE)

args_list = argsReshape(default_args, args, num_args, log_args)

print(args_list)


print("** Loading packages and data... **")
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(Matrix, quietly = T))
suppressPackageStartupMessages(library(seqminer, quietly = T))
suppressPackageStartupMessages(library(SurvSTAAR, quietly = T))


print(sessionInfo())


rap_load_as(args_list$objnull.file, "objNull")
print(paste0("Load objnull data successfully! ", args_list$objnull.file))



# individual analysis -----------------------------------------------------


output_name = paste0(args_list$output.file,".rda")

bim_data = fread(paste0(args_list$plink.file, ".bim"))
length_allSet = nrow(bim_data)

begin = Sys.time()
print(paste0("Begin running the program at ", begin))

if (is.null(args_list$set) & is.null(args_list$Nsets) & is.null(args_list$start_loc) & is.null(args_list$end_loc)) {
  ## analyze all variants in plink file in one job

  print(paste0("In this job, there are ", length_allSet, " single variants will be tested"))

  individual_results = IndividualTestPlink(objNull = objNull,
                                           genofile = args_list$plink.file,
                                           chr = args_list$chr,
                                           sampleCol =  args_list$sampleCol,
                                           use_SPA = args_list$use_SPA,
                                           SPA_filter = args_list$SPA_filter,
                                           SPA_filter_cutoff = args_list$SPA_filter_cutoff,
                                           geno_missing_cutoff = args_list$geno_missing_cutoff,
                                           geno_missing_imputation = args_list$geno_missing_imputation,
                                           min_mac_cutoff = args_list$min_mac_cutoff,
                                           min_maf_cutoff = args_list$min_maf_cutoff,
                                           chunk_size = args_list$chunk_size,
                                           verbose = args_list$verbose)

} else if (is.null(args_list$set) & is.null(args_list$Nsets) & !is.null(args_list$start_loc) & !is.null(args_list$end_loc)) {
  ## analyze all user defined variants within start_loc and end_loc

  individual_results = IndividualTestPlink(objNull = objNull,
                                           start_loc = args_list$start_loc,
                                           end_loc = args_list$end_loc,
                                           genofile = args_list$plink.file,
                                           chr = args_list$chr,
                                           sampleCol =  args_list$sampleCol,
                                           use_SPA = args_list$use_SPA,
                                           SPA_filter = args_list$SPA_filter,
                                           SPA_filter_cutoff = args_list$SPA_filter_cutoff,
                                           geno_missing_cutoff = args_list$geno_missing_cutoff,
                                           geno_missing_imputation = args_list$geno_missing_imputation,
                                           min_mac_cutoff = args_list$min_mac_cutoff,
                                           min_maf_cutoff = args_list$min_maf_cutoff,
                                           chunk_size = args_list$chunk_size,
                                           verbose = args_list$verbose)

} else if (!is.null(args_list$set) & !is.null(args_list$Nsets) & is.null(args_list$start_loc) & is.null(args_list$end_loc)) {
  ## separate all variants into "Nsets" jobs, the start_loc and end_loc are defined automatically

  Nsets = args_list$Nsets
  set = args_list$set

  if (Nsets != 1) {
    print(paste0("All the masks will be devided into ", Nsets, " groups"))

    if (set == Nsets) {
      numSet = length_allSet - ceiling(length_allSet / Nsets) * (Nsets - 1)
      bim_data_Set = bim_data[(length_allSet - numSet + 1):length_allSet, ]
      args_list$start_loc = min(bim_data_Set[, 4])
      args_list$end_loc = max(bim_data_Set[, 4])
    } else {
      numSet = ceiling(length_allSet / Nsets)
      bim_data_Set = bim_data[(numSet * (set - 1) + 1):(numSet * set), ]
      args_list$start_loc = min(bim_data_Set[, 4])
      args_list$end_loc = max(bim_data_Set[, 4])
    }

    print(paste0("In sub-job ", set , ", there are ", numSet, " masks will be tested"))

    individual_results = IndividualTestPlink(objNull = objNull,
                                             start_loc = args_list$start_loc,
                                             end_loc = args_list$end_loc,
                                             genofile = args_list$plink.file,
                                             chr = args_list$chr,
                                             sampleCol =  args_list$sampleCol,
                                             use_SPA = args_list$use_SPA,
                                             SPA_filter = args_list$SPA_filter,
                                             SPA_filter_cutoff = args_list$SPA_filter_cutoff,
                                             geno_missing_cutoff = args_list$geno_missing_cutoff,
                                             geno_missing_imputation = args_list$geno_missing_imputation,
                                             min_mac_cutoff = args_list$min_mac_cutoff,
                                             min_maf_cutoff = args_list$min_maf_cutoff,
                                             chunk_size = args_list$chunk_size,
                                             verbose = args_list$verbose)


  } else if (Nsets == 1) {

    print(paste0("In this job, there are ", length_allSet, " single variants will be tested"))

    individual_results = IndividualTestPlink(objNull = objNull,
                                             genofile = args_list$plink.file,
                                             chr = args_list$chr,
                                             sampleCol =  args_list$sampleCol,
                                             use_SPA = args_list$use_SPA,
                                             SPA_filter = args_list$SPA_filter,
                                             SPA_filter_cutoff = args_list$SPA_filter_cutoff,
                                             geno_missing_cutoff = args_list$geno_missing_cutoff,
                                             geno_missing_imputation = args_list$geno_missing_imputation,
                                             min_mac_cutoff = args_list$min_mac_cutoff,
                                             min_maf_cutoff = args_list$min_maf_cutoff,
                                             chunk_size = args_list$chunk_size,
                                             verbose = args_list$verbose)

  }

}


end = Sys.time()
print(paste0("Finish analyzing this job at ", end))
elapse = difftime(end, begin, units = "hours")

print(paste0("The results has been saved at ", output_name))
