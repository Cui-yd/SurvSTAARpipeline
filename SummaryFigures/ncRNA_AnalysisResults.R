##########################################################
# Summarize and visualize ncRNA masks analysis results
# Yidan Cui
# Initiate date: 2025/02/28
# Current date: 2025/02/28
##########################################################

setwd("/PATH/TO/RESULTS")

ncRNA_all = NULL

file = list.files(path = paste0("./"), pattern = ".rda")
filename = paste0("./", file)
pb <- txtProgressBar(min = 1, max = length(filename), style = 3)


for (i in 1:length(filename)) {
  rap_load_as(filename[i], "ncRNA_temp")
  ncRNA_all = c(ncRNA_all, ncRNA_temp)
  setTxtProgressBar(pb, i)
}


survstaar_ncRNA = NULL
pb <- txtProgressBar(min = 1, max = length(ncRNA_all), style = 3)

for (j in 1:length(ncRNA_all)) {

  if (inherits(ncRNA_all[[j]], "try-error")) {
    survstaar_ncRNA_temp = rep(NA, 3)
  } else {
    survstaar_ncRNA_temp = as.data.frame(t(ncRNA_all[[j]]$SurvSTAAR_O))
    survstaar_ncRNA_temp=cbind(ncRNA_all[[j]]$gene_info, survstaar_ncRNA_temp)
  }

  survstaar_ncRNA = rbind(survstaar_ncRNA, survstaar_ncRNA_temp)

  setTxtProgressBar(pb, j)
}

table(survstaar_ncRNA$Chr)
table(ncRNA_info$Chr)

no_results = setdiff(ncRNA_info$ncRNA, survstaar_ncRNA$ncRNA)
ncRNA_info[which(ncRNA_info$ncRNA %in% no_results),]

survstaar_ncRNA = survstaar_ncRNA[na.omit(match(ncRNA_info$ncRNA, survstaar_ncRNA$ncRNA)),]
ncRNA_info$mid_pos = (ncRNA_info$Start_POS + ncRNA_info$End_POS) / 2
survstaar_ncRNA$pos = ncRNA_info$mid_pos[match(survstaar_ncRNA$ncRNA, ncRNA_info$ncRNA)]

save(survstaar_ncRNA, file = "./survstaar_ncRNA.rda")


## Manhattan plot
colnames(survstaar_ncRNA)[3] = c("ncRNA")

png(filename = "./manhattan_ncRNA.png", width = 13, height = 8, units = "in", res = 300)
SurvSTAAR_Manhattan(result_data = survstaar_ncRNA,
                    name = "ncRNA", chr = "Chr", pos = "pos",
                    col = c("#D68800", "#12008B"), max_y = 20,
                    pval = c("ncRNA"), pch = 1, 
                    annotateP = 2.5e-6, genoline = 2.5e-6, legned_title = "")

## no label
# SurvSTAAR_Manhattan(result_data = survstaar_ncRNA,
#                     name = "ncRNA", chr = "Chr", pos = "pos",
#                     col = c("#D68800", "#12008B"), max_y = 20,
#                     pval = c("ncRNA"), pch = 1, 
#                     genoline = 2.5e-6, legned_title = "")
dev.off()


## QQplot
pval_result = survstaar_ncRNA[,3,F]

png(filename = "qqplot_ncRNA.png", width = 7, height = 7, units = "in", res = 300)
SurvSTAAR_QQplot(pval_result = pval_result, pch = 1, legned_title = FALSE,
                 legend_lable = c("ncRNA"))
dev.off()
