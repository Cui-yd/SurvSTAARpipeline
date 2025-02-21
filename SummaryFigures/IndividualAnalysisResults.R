##########################################################
# Summarize and visualize individual analysis results
# Yidan Cui
# Initiate date: 2025/02/21
# Current date: 2025/02/21
##########################################################

setwd("/PATH/TO/RESULTS")

survstaar_individual = NULL

file = list.files(path = ".", pattern = "rda")
filename = paste0("./", file)
pb <- txtProgressBar(min = 1, max = length(filename), style = 3)

for (i in 1:length(filename)) {
  rap_load_as(filename[i], "survstaar_individual_temp")
  survstaar_individual = rbind(survstaar_individual, survstaar_individual_temp)
  setTxtProgressBar(pb, i)
}

print(dim(survstaar_individual))

## Manhattan plot
SurvSTAAR_Manhattan(result_data = survstaar_individual,
                    name = "rsID", chr = "CHR", pos = "POS",
                    col = c("#D68800", "#12008B"), max_y = 20, pval = "Pvalue",
                    pch = 1, genoline = 5e-8, legned_title = "")


## QQplot
pval_result = survstaar_individual$Pvalue_SPA
SurvSTAAR_QQplot(pval_result = pval_result, pch = 1, legned_title = FALSE, legend_lable = "")

