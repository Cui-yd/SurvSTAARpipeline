##########################################################
# Analyze noncoding functional categories using SurvSTAAR
# on High Performance Computing Cluster
# Yidan Cui
# Initiate date: 2025/02/20
# Current date: 2025/04/10
##########################################################

rm(list = ls())
gc()

## Load packages

print("** Loading packages and data... **")
suppressPackageStartupMessages(library(gdsfmt, quietly = T))
suppressPackageStartupMessages(library(SeqArray, quietly = T))
suppressPackageStartupMessages(library(SeqVarTools, quietly = T))
suppressPackageStartupMessages(library(Matrix, quietly = T))
suppressPackageStartupMessages(library(stats, quietly = T))
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(CompQuadForm, quietly = T))
suppressPackageStartupMessages(library(SurvSTAAR, quietly = T))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = T))

print(sessionInfo())


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
  categories = "all",
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
if (length(strsplit(args_list$Annotation_name, ",")) == 1) {
  args_list$Annotation_name = strsplit(args_list$Annotation_name, ",")[[1]]
}

print(args_list)


genofile = seqOpen(args_list$agds.file)
print(paste0("Open aGDS file successfully! ", basename(genofile$filename)))

rap_load_as(args_list$objnull.file, "objNull")
print(paste0("Load objnull data successfully! ", args_list$objnull.file))

rap_load_as(args_list$Annotation.name.catalog.file, "Annotation_name_catalog")
print(paste0("Load annotation catalog data successfully! ", args_list$Annotation.name.catalog.file))


genes_info_chr = genes_info[genes_info$Chr == args_list$chr, ]
print(paste0("Number of masks in chr", args_list$chr, " is: ", nrow(genes_info_chr)))


Nsets = args_list$Nsets
set = args_list$set
length_allSet = nrow(genes_info_chr)

if (Nsets != 1) {
  print(paste0("All the masks will be devided into ", args_list$Nsets, " groups"))

  if (set == Nsets) {
    numSet = length_allSet - ceiling(length_allSet / Nsets) * (Nsets - 1)
    idSet = genes_info_chr$Gene[(length_allSet - numSet + 1):length_allSet]
  } else {
    numSet = ceiling(length_allSet / Nsets)
    idSet = genes_info_chr$Gene[(numSet * (set - 1) + 1):(numSet * set)]
  }

  print(paste0("In sub-job ", set , ", there are ", numSet, " masks will be tested"))

} else if (Nsets == 1) {

  numSet = length_allSet
  print(paste0("In this job, there are ", numSet, " masks will be tested"))

}



# preload information -----------------------------------------------------


### Promoter_CAGE ######

varid <- seqGetData(genofile, "variant.id")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promGobj <- suppressMessages(promoters(genes(txdb), upstream = 3000, downstream = 3000))

#Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
CAGEBvt <- CAGEAnno!=""
CAGEidx <- which(CAGEBvt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[CAGEidx])
seqSetFilter(genofile,promGobj,intersect=TRUE)
CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
##obtain variants info
CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
CAGEvref <- as.character(seqGetData(genofile,"$ref"))
CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant")
{
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV")
{
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel")
{
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.PromCAGE <- variant.id[SNVlist]

dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

seqResetFilter(genofile)

rm(dfPromCAGEVarGene)
gc()


### Promoter_DHS ######

varid <- seqGetData(genofile, "variant.id")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

# Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRsBvt <- rOCRsAnno!=""
rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[rOCRsidx])

seqSetFilter(genofile,promGobj,intersect=TRUE)
rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
## obtain variants info
rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant")
{
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV")
{
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel")
{
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.PromrOCRs <- variant.id[SNVlist]

dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

seqResetFilter(genofile)

rm(dfPromrOCRsVarGene)
gc()


### Enhancer_CAGE ######

varid <- seqGetData(genofile, "variant.id")

#Now extract the GeneHancer with CAGE Signal Overlay
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
genehancer <- genehancerAnno!=""

CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
CAGE <- CAGEAnno!=""
CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

# variants that covered by whole GeneHancer without CAGE overlap.
genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
enhancervpos <- as.numeric(seqGetData(genofile,"position"))
enhancervref <- as.character(seqGetData(genofile,"$ref"))
enhancervalt <- as.character(seqGetData(genofile,"$alt"))
dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant")
{
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV")
{
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel")
{
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.HancerCAGE <- variant.id[SNVlist]

dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)

seqResetFilter(genofile)

rm(dfHancerCAGEVarGene)
gc()


### Enhancer_DHS ######

varid <- seqGetData(genofile, "variant.id")

#Now extract the GeneHancer with rOCRs Signal Overlay
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
genehancer <- genehancerAnno!=""

rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRs <- rOCRsAnno!=""
rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
# variants that covered by whole GeneHancer without rOCRs overlap.

genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
enhancervpos <- as.numeric(seqGetData(genofile,"position"))
enhancervref <- as.character(seqGetData(genofile,"$ref"))
enhancervalt <- as.character(seqGetData(genofile,"$alt"))
dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant")
{
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV")
{
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel")
{
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.HancerrOCRs <- variant.id[SNVlist]

dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

seqResetFilter(genofile)

rm(dfHancerrOCRsVarGene)
gc()


print("** All the data and functions required are settled down! **")



# gene-centric noncoding analysis --------------------------------------------

output_name = paste0(args_list$output.file,".rda")

geneNonCoding_results = list()
time_elapse = c()
begin = Sys.time()

print(paste0("Begin running the program at ", begin))

for (j in 1:length(idSet)) {

  begin_j = Sys.time()

  gene_name = idSet[j]
  print(paste0(j, "/", length(idSet), "      ", gene_name))

  geneNonCoding_results[[j]] = try(geneCentricNonCoding_preload(gene_name = gene_name,
                                                                genofile = genofile,
                                                                objNull = objNull,
                                                                genes_info = genes_info_chr,
                                                                variant_type = args_list$variant_type,
                                                                categories = args_list$categories,
                                                                dfPromCAGEVarGene.SNV = dfPromCAGEVarGene.SNV,
                                                                variant.id.SNV.PromCAGE = variant.id.SNV.PromCAGE,
                                                                dfPromrOCRsVarGene.SNV = dfPromrOCRsVarGene.SNV,
                                                                variant.id.SNV.PromrOCRs = variant.id.SNV.PromrOCRs,
                                                                dfHancerCAGEVarGene.SNV = dfHancerCAGEVarGene.SNV,
                                                                variant.id.SNV.HancerCAGE = variant.id.SNV.HancerCAGE,
                                                                dfHancerrOCRsVarGene.SNV = dfHancerrOCRsVarGene.SNV,
                                                                variant.id.SNV.HancerrOCRs = variant.id.SNV.HancerrOCRs,
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

  save(geneNonCoding_results, file = output_name)

}


end = Sys.time()
print(paste0("Finish analyzing this job at ", end))
elapse = difftime(end, begin, units = "hours")
print(paste0("Time elapse: ", round(elapse, 2), " hours"))
print(paste0("Average time consuming per mask is ", round(mean(time_elapse), 2), " mins"))

seqClose(genofile)
gc()

print(paste0("The results has been saved at ", output_name))
