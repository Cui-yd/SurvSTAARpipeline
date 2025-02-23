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

R/Bioconductor package **SeqArray** provides functions to convert the genotype data (in VCF/BCF/PLINK BED/SNPRelate format) to SeqArray GDS format. For more details on usage, please see the R/Bioconductor package <a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">**SeqArray**</a> [<a href="https://www.bioconductor.org/packages/devel/bioc/manuals/SeqArray/man/SeqArray.pdf">manual</a>].

R package **gds2bgen** provides functions to convert the genotype data (in BGEN format) to SeqArray GDS format. For more details on usage, please see the R package <a href="https://github.com/zhengxwen/gds2bgen">**gds2bgen**</a>. An example for the `seqBGEN2GDS` function in the gds2bgen package can be found <a href="https://github.com/zhengxwen/gds2bgen#examples">**here**</a> (**Credit: Xiuwen Zheng**).

**FAVORannotator** (CSV version 1.0.0) depends on the **xsv software** and the **FAVOR database** in CSV format. Please install the <a href="https://github.com/BurntSushi/xsv">**xsv software**</a> and download the **FAVOR essential database CSV files** from <a href="http://favor.genohub.org">**FAVOR website**</a> (under the "FAVORannotator" tab's top panel, 31.2 GB for chr1 CSV) or <a href="https://doi.org/10.7910/DVN/1VGTJI">**Harvard Dataverse**</a> before using **FAVORannotator** (CSV version 1.0.0).

For step-by-step instructions, please refer to <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">**STAARpipeline-Tutorial**</a> and <a href="http://favor.genohub.org">**FAVOR website**</a> for the relevant steps.


### Step 0: Calculate polygenic effects

In this step, we highly recommend users to use genotype array data (originally in PLINK format) to calculate polygenic effects. The corresponding Bash scripts can be found in the `script_step0` folder.

#### Prepare the genotype data
- Use <a href="https://www.cog-genomics.org/plink/">**PLINK1.9**</a> the 22 chromosome PLINK files into a single file using `script_step0/merge.sh`.
- Use <a href="https://www.cog-genomics.org/plink/2.0/">**PLINK2**</a> to perform pruning on the merged file with the specific coefficients recommended by <a href="https://www.nature.com/articles/s41588-021-00870-7"> **REGENIE (Mbatchou, J. et al.)** </a>, using `script_step0/pruen.sh`.
- Use <a href="https://www.cog-genomics.org/plink/2.0/">**PLINK2**</a> to extract the selected variants into a new PLINK file using `script_step0/`.
#### Prepare the phenotype
- Please refer to <a href="https://github.com/Cui-yd/ukbSurvPhe">**ukbSurvPhe**</a> for more information on setting up your time-to-event phenotype.
#### Inpute:
- Pruned genotype array data in PLINK format.
- Phenotype data containing only survival status (with header: IID FID status).
- Covariate data (with header: IID FID cov1 cov2 ...).
#### Submit the job
In this step, we calculate polygenic effects using a whole-genome regression model with REGENIE. Users can refer to `script_step0/submit.sh` for submitting the job through REGENIE. For more details, including installation and additional job submission parameter settings, please refer to the <a href="https://rgcgithub.github.io/regenie/overview/#step-1-whole-genome-model">**REGENIE Step 1**</a> documentation.
#### Output:
- `regenie_result_1.loco`
- `regenie_result.log`

### Step 1: Fit Cox proportional hazards null model

#### Prepare the phenotype for null model
Users can combine the phenotype, covariates, and polygenic effects into a new data table using script `script_step1/combine_phenotype_polygenic.R`.
##### Inpute:
- Phenotype data in .txt format (with header: IID FID time status).
- Output from step 0 `regenie_result_1.loco`.
- Covariate data (with header: IID FID cov1 cov2 ...).
##### Output:
- `/phenotype_all.txt` (with header: IID FID time status #Chr1:22... #all_covariates...)
  
#### Fit null model
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
