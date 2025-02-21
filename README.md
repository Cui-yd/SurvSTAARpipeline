# SurvSTAARpipeline
The analytical pipeline and implementation for the UK Biobank Research Analysis Platform (RAP) of SurvSTAAR


## Description
SurvSTAAR is an R package for performing variant-Set Test for Association using Annotation infoRmation (STAAR) procedure in whole-genome sequencing (WGS) studies in time-to-event traits. 
SurvSTAAR is a general framework that incorporates both qualitative functional categories and quantitative complementary functional annotations using an omnibus test SurvSTAAR-O. 
SurvSTAAR accounts for population structure and sample relatedness, and addresses challenges posed by heavily censored phenotypes and low-frequency variants.



## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 4.0.0)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library 
(or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.


## Dependencies

SurvSTAAR imports R packages 
<a href="https://cran.r-project.org/web/packages/survival/index.html">survival</a>,
<a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>, 
<a href="https://cran.r-project.org/web/packages/data.table/index.html">data.table</a>,
<a href="https://cran.r-project.org/web/packages/CompQuadForm/index.html">CompQuadForm</a>,
<a href="https://cran.r-project.org/web/packages/dplyr/index.html">dplyr</a>,
<a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">SeqArray</a>,
<a href="https://bioconductor.org/packages/release/bioc/html/SeqVarTools.html">SeqVarTools</a>,
<a href="https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html">GenomicFeatures</a>,
<a href="https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html">TxDb.Hsapiens.UCSC.hg38.knownGene</a>,
 
These dependencies should be installed before using SurvSTAAR and this pipeline.


## Installation
```
library(devtools)
devtools::install_github("Cui-yd/SurvSTAAR",ref="main")
```

## Usage

### Step 0: Generate annotated GDS (aGDS) file using FAVORannotator

Please refer to <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">STAARpipeline-Tutorial</a> for the relevant steps.


### Step 0: Calculate polygenic effects

Please refer to <a href="https://rgcgithub.github.io/regenie/overview/#step-1-whole-genome-model">REGENIE-step1</a> for the relevant steps.

### Step 1: Fit Cox proportional hazards null model

`fitNullModel_onCluster.R` and `fitNullModel_onRAP.R` scripts designed for the HPC and RAP platforms, respectively. 
You can use the `generate_NullModel_command.sh` script to generate the corresponding execution command. 
If you are using the RAP platform, please submit the job via Swiss Army Knife.


### Step 2a: Individual analysis for common varaints

`IndividualAnalysisGDS_onCluster.R` and `IndividualAnalysisGDS_onRAP.R` scripts designed for the aGDS format files to perform individual analysis.
`IndividualAnalysisPlink_onCluster.R` and `IndividualAnalysisPlink_onRAP.R` scripts designed for the plink format files to perform individual analysis.
You can use the `generate_IndividualGDS_command.sh` and `generate_IndividualPlink_command.sh` script to generate the corresponding execution command. 
If you are using the RAP platform, please submit the job via Swiss Army Knife.



### Step 2b: Gene-based test for rare variants

#### Gene-centric coding analysis
`GeneCentricCoding_onCluster.R` and `GeneCentricCoding_onRAP.R` scripts designed for the aGDS format files to perform gene-centric coding analysis.
You can use the `generate_Coding_command.sh` script to generate the corresponding execution command. 
If you are using the RAP platform, please submit the job via Swiss Army Knife.

#### Gene-centric noncoding analysis
`GeneCentricNonCoding_onRAP.R` and `GeneCentricNonCoding_onCluster.R` scripts designed for the aGDS format files to perform gene-centric noncoding analysis.
You can use the `generate_NonCoding_command.sh` script to generate the corresponding execution command. 
If you are using the RAP platform, please submit the job via Swiss Army Knife.

`ncRNA_onCluster.R` and `ncRNA_onRAP.R` scripts designed for the aGDS format files to perform gene-centric noncoding analysis for ncRNA genes across the genome.
You can use the `generate_ncRNA_command.sh` script to generate the corresponding execution command. 
If you are using the RAP platform, please submit the job via Swiss Army Knife.


### Data processing and visualization
`CodingAnalysisResults.R`, `IndividualAnalysisResults.R`, and `NonCodingAnalysisResults.R` are designed for organizing the corresponding analysis results and generating Manhattan and QQ plots.



## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
