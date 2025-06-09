##########################################################
# Fit Cox proportional hazards null model using SurvSTAAR
# on UK Biobank Research Analysis Platform
# Yidan Cui
# Initiate date: 2025/02/20
# Current date: 2025/04/10
##########################################################

rm(list = ls())
gc()

## Load packages

suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(survival, quietly = T))
suppressPackageStartupMessages(library(Matrix, quietly = T))
suppressPackageStartupMessages(library(SurvSTAAR, quietly = T))

print(sessionInfo())

# set parameters ----------------------------------------------------------

default_args = list(
  chr = 666,
  LOCO = TRUE,
  genofile = NULL,
  genofile_class = NULL,
  phenofile = "",
  statusCol = "status",
  timeCol = "time",
  covCol = c("sex", "birthyr", "PC1", "PC2", "PC3", "PC4", "PC5",
             "PC6", "PC7", "PC8", "PC9", "PC10"),
  sampleCol = "IID",
  output = "output_default",
  use_SPA = TRUE,
  verbose = FALSE
)

num_args = "chr"
log_args = c("LOCO", "use_SPA", "verbose")

args <- commandArgs(TRUE)

args_list = argsReshape(default_args, args, num_args, log_args)
if (length(strsplit(args_list$covCol, ",")) == 1) {
  args_list$covCol = strsplit(args_list$covCol, ",")[[1]]
}

print(args_list)


# fit null model -----------------------------------------------------------

if (!is.null(args_list$genofile)) {
  if (args_list$genofile_class == "agds") {
    genofile = seqOpen(args_list$genofile)
  }
}


PRSCol = paste0("Chr", args_list$chr)
output_name = paste0(args_list$output, ".rda")

begin = Sys.time()

objNull = NullModel(genofile = args_list$genofile, phenofile = args_list$phenofile,
                    LOCO = args_list$LOCO, chr = args_list$chr,
                    statusCol = args_list$statusCol, timeCol = args_list$timeCol,
                    sampleCol = args_list$sampleCol, covCol = args_list$covCol, PRSCol = PRSCol,
                    use_SPA = args_list$use_SPA, range=c(-100,100), length.out = 1e4,
                    verbose = args_list$verbose)

print(paste0("Elapsed time is: ", round(difftime(Sys.time(), begin, units = "secs"), 2), " secs"))

save(objNull, file = output_name)

if (!is.null(args_list$genofile)) {
  if (args_list$genofile_class == "agds") {
    seqClose(genofile)
    gc()
  }
}


print(paste0("The results has been saved at ", output_name))


system(paste0("dx upload ", output_name))
