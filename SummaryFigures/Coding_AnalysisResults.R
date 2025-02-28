##########################################################
# Summarize and visualize coding masks analysis results
# Yidan Cui
# Initiate date: 2025/02/21
# Current date: 2025/02/28
##########################################################

setwd("/PATH/TO/RESULTS")

geneCoding_all = NULL

file = list.files(path = paste0("./"), pattern = ".rda")
filename = paste0("./", file)
pb <- txtProgressBar(min = 1, max = length(filename), style = 3)


for (i in 1:length(filename)) {
  rap_load_as(filename[i], "geneCoding_temp")
  geneCoding_all = c(geneCoding_all, geneCoding_temp)
  setTxtProgressBar(pb, i)
}


survstaar_coding = NULL
pb <- txtProgressBar(min = 1, max = length(geneCoding_all), style = 3)

for (j in 1:length(geneCoding_all)) {

  if (inherits(geneCoding_all[[j]], "try-error")) {
    survstaar_coding_temp = rep(NA, 7)
  } else {
    survstaar_coding_temp = as.data.frame(t(geneCoding_all[[j]]$SurvSTAAR_O_all))
    survstaar_coding_temp=cbind(geneCoding_all[[j]]$gene_info, survstaar_coding_temp)
  }

  survstaar_coding = rbind(survstaar_coding, survstaar_coding_temp)

  setTxtProgressBar(pb, j)
}

table(survstaar_coding$Chr)
table(genes_info$Chr)

no_results = setdiff(genes_info$Gene, survstaar_coding$Gene)
genes_info[which(genes_info$Gene %in% no_results),]

survstaar_coding = survstaar_coding[na.omit(match(genes_info$Gene, survstaar_coding$Gene)),]
genes_info$mid_pos = (genes_info$Start_POS + genes_info$End_POS) / 2
survstaar_coding$pos = genes_info$mid_pos[match(survstaar_coding$Gene, genes_info$Gene)]

save(survstaar_coding, file = "./survstaar_coding.rda")


## Manhattan plot
colnames(survstaar_coding)[3:7] = c("pLoF", "pLoF+D", "Missense", "Disruptive Missense", "Synonymous")

png(filename = "./manhattan_coding.png", width = 13, height = 8, units = "in", res = 300)
SurvSTAAR_Manhattan(result_data = survstaar_coding,
                    name = "Gene", chr = "Chr", pos = "pos",
                    col = c("#D68800", "#12008B"), max_y = 20,
                    pval = c("pLoF", "pLoF+D", "Missense", "Disruptive Missense", "Synonymous"),
                    pch = 0:4, annotateP = 5e-7, genoline = 5e-7, legned_title = "Functional Categories")

## no label
# SurvSTAAR_Manhattan(result_data = survstaar_coding,
#                     name = "Gene", chr = "Chr", pos = "pos",
#                     col = c("#D68800", "#12008B"), max_y = 20,
#                     pval = c("pLoF", "pLoF+D", "Missense", "Disruptive Missense", "Synonymous"),
#                     pch = 0:4, genoline = 5e-7, legned_title = "Functional Categories")
dev.off()


## QQplot
pval_result = survstaar_coding[,3:7]

png(filename = "qqplot_coding.png", width = 7, height = 7, units = "in", res = 300)
SurvSTAAR_QQplot(pval_result = pval_result, pch = 0:4, legned_title = FALSE,
                 legend_lable = c("pLoF", "pLoF+D", "Missense", "Disruptive Missense", "Synonymous"))
dev.off()
