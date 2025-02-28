##########################################################
# Summarize and visualize noncoding masks analysis results
# Yidan Cui
# Initiate date: 2025/02/21
# Current date: 2025/02/28
##########################################################

setwd("/PATH/TO/RESULTS")

geneNonCoding_all = NULL

file = list.files(path = paste0("./"), pattern = ".rda")
filename = paste0("./", file)
pb <- txtProgressBar(min = 1, max = length(filename), style = 3)


for (i in 1:length(filename)) {
  rap_load_as(filename[i], "geneNonCoding_temp")
  geneNonCoding_all = c(geneNonCoding_all, geneNonCoding_temp)
  setTxtProgressBar(pb, i)
}


survstaar_noncoding = NULL
pb <- txtProgressBar(min = 1, max = length(geneNonCoding_all), style = 3)

for (j in 1:length(geneNonCoding_all)) {

  if (inherits(geneNonCoding_all[[j]], "try-error")) {
    survstaar_noncoding_temp = rep(NA, 9)
  } else {
    survstaar_noncoding_temp = as.data.frame(t(geneNonCoding_all[[j]]$SurvSTAAR_O_all))
    survstaar_noncoding_temp=cbind(geneNonCoding_all[[j]]$gene_info, survstaar_noncoding_temp)
  }

  survstaar_noncoding = rbind(survstaar_noncoding, survstaar_noncoding_temp)

  setTxtProgressBar(pb, j)
}

table(survstaar_noncoding$Chr)
table(genes_info$Chr)

no_results = setdiff(genes_info$Gene, survstaar_noncoding$Gene)
genes_info[which(genes_info$Gene %in% no_results),]

survstaar_noncoding = survstaar_noncoding[na.omit(match(genes_info$Gene, survstaar_noncoding$Gene)),]
genes_info$mid_pos = (genes_info$Start_POS + genes_info$End_POS) / 2
survstaar_noncoding$pos = genes_info$mid_pos[match(survstaar_noncoding$Gene, genes_info$Gene)]

save(survstaar_noncoding, file = "./survstaar_noncoding.rda")


## Manhattan plot
colnames(survstaar_noncoding)[3:9] = c("Upstream","Downstream","UTR","Promoter_CAGE","Promoter_DHS","Enhancer_CAGE","Enhancer_DHS")

png(filename = "./manhattan_noncoding.png", width = 13, height = 8, units = "in", res = 300)
SurvSTAAR_Manhattan(result_data = survstaar_noncoding,
                    name = "Gene", chr = "Chr", pos = "pos",
                    col = c("#D68800", "#12008B"), max_y = 20,
                    pval = c("Upstream","Downstream","UTR","Promoter_CAGE","Promoter_DHS","Enhancer_CAGE","Enhancer_DHS"),
                    pch = 0:6, annotateP = 3.57e-7, genoline = 3.57e-7, legned_title = "Functional Categories")

## no label
# SurvSTAAR_Manhattan(result_data = survstaar_noncoding,
#                     name = "Gene", chr = "Chr", pos = "pos",
#                     col = c("#D68800", "#12008B"), max_y = 20,
#                     pval = c("Upstream","Downstream","UTR","Promoter_CAGE","Promoter_DHS","Enhancer_CAGE","Enhancer_DHS"),
#                     pch = 0:6, genoline = 3.57e-7, legned_title = "Functional Categories")
dev.off()


## QQplot
pval_result = survstaar_noncoding[,3:9]

png(filename = "qqplot_noncoding.png", width = 7, height = 7, units = "in", res = 300)
SurvSTAAR_QQplot(pval_result = pval_result, pch = 0:6, legned_title = FALSE,
                 legend_lable = c("Upstream","Downstream","UTR","Promoter_CAGE","Promoter_DHS","Enhancer_CAGE","Enhancer_DHS"))
dev.off()
