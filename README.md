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
<a href="https://cran.r-project.org/web/packages/seqminer/index.html">seqminer</a>,
<a href="https://cran.r-project.org/web/packages/dplyr/index.html">dplyr</a>,
<a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">SeqArray</a>,
<a href="https://bioconductor.org/packages/release/bioc/html/SeqVarTools.html">SeqVarTools</a>,
<a href="https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html">GenomicFeatures</a>,
<a href="https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html">TxDb.Hsapiens.UCSC.hg38.knownGene</a>.
 
These dependencies should be installed before using SurvSTAAR and this pipeline.


## Installation
```
install.packages(c("survival", "Matrix", "data.table", "CompQuadForm", "seqminer", "dplyr"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SeqArray", "SeqVarTools", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene"))

devtools::install_github("Cui-yd/SurvSTAAR", ref = "main")
```

## Usage

### Step 0.1: Generate annotated GDS (aGDS) file using FAVORannotator

(**Skip this step if you are only analyzing common variants at the individual level.**)

R/Bioconductor package **SeqArray** provides functions to convert the genotype data (in VCF/BCF/PLINK BED/SNPRelate format) to SeqArray GDS format. For more details on usage, please see the R/Bioconductor package <a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">**SeqArray**</a> [<a href="https://www.bioconductor.org/packages/devel/bioc/manuals/SeqArray/man/SeqArray.pdf">manual</a>].

R package **gds2bgen** provides functions to convert the genotype data (in BGEN format) to SeqArray GDS format. For more details on usage, please see the R package <a href="https://github.com/zhengxwen/gds2bgen">**gds2bgen**</a>. An example for the `seqBGEN2GDS` function in the gds2bgen package can be found <a href="https://github.com/zhengxwen/gds2bgen#examples">**here**</a> (**Credit: Xiuwen Zheng**).

**FAVORannotator** (CSV version 1.0.0) depends on the **xsv software** and the **FAVOR database** in CSV format. Please install the <a href="https://github.com/BurntSushi/xsv">**xsv software**</a> and download the **FAVOR essential database CSV files** from <a href="http://favor.genohub.org">**FAVOR website**</a> (under the "FAVORannotator" tab's top panel, 31.2 GB for chr1 CSV) or <a href="https://doi.org/10.7910/DVN/1VGTJI">**Harvard Dataverse**</a> before using **FAVORannotator** (CSV version 1.0.0).

For step-by-step instructions, please refer to <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">**STAARpipeline-Tutorial**</a> and <a href="http://favor.genohub.org">**FAVOR website**</a> for the relevant steps.


### Step 0.2: Calculate polygenic effects

In this step, we highly recommend users to use genotype array data (originally in PLINK format) to calculate polygenic effects. The corresponding Bash scripts can be found in the **`script_step0`** folder.

#### Step 0.2.1 Prepare the genotype array data
- Use <a href="https://www.cog-genomics.org/plink/">**PLINK1.9**</a> merge the 22 chromosome PLINK files into a single file using <a href="script_step0/merge.sh">**script_step0/merge.sh**</a>.
- Use <a href="https://www.cog-genomics.org/plink/2.0/">**PLINK2**</a> to perform pruning on the merged file with the specific coefficients recommended by <a href="https://www.nature.com/articles/s41588-021-00870-7"> **REGENIE (Mbatchou, J. et al.)** </a>, using <a href="script_step0/pruen.sh">**script_step0/pruen.sh**</a>.
- Use <a href="https://www.cog-genomics.org/plink/2.0/">**PLINK2**</a> to extract the selected variants into a new PLINK file using <a href="script_step0/extract.sh">**script_step0/extract.sh**</a>.
  
#### Step 0.2.2 Prepare the phenotype
- Please refer to <a href="https://github.com/Cui-yd/ukbSurvPhe">**ukbSurvPhe**</a> for more information on setting up your time-to-event phenotype.
##### Input:
- Pruned genotype array data in PLINK format.
- Phenotype data containing only survival status (with header: IID FID status).
- Covariate data (with header: IID FID cov1 cov2 ...).
##### Submit the job
In this step, we calculate polygenic effects using a whole-genome regression model with REGENIE. Users can refer to <a href="script_step0/submit.sh">**script_step0/submit.sh**</a> for submitting the job through REGENIE. For more details, including installation and additional job submission parameter settings, please refer to the <a href="https://rgcgithub.github.io/regenie/overview/#step-1-whole-genome-model">**REGENIE Step 1**</a> documentation. (For reference, this step takes approximately 2 hours for 400,000 genotype variants.)
##### Output:
- `regenie_result_1.loco`
- `regenie_result.log`


### Step 1: Fit Cox proportional hazards null model

#### Step 1.1 Prepare the phenotype for null model
Users can combine the phenotype, covariates, and polygenic effects into a new data table using script <a href="script_step1/combine_phenotype_polygenic.R">**script_step1/combine_phenotype_polygenic.R**</a>.
##### Input:
- Phenotype data in .txt format (with header: IID FID time status).
- Output from step 0 `regenie_result_1.loco`.
- Covariate data (with header: IID FID cov1 cov2 ...).
##### Output:
- `phenotype_all.txt` (with header: IID FID time status #Chr1:22... #all_covariates...)
  
#### Step 1.2 Fit null model
<a href="HPC_script/fitNullModel_onCluster.R">**HPC_script/fitNullModel_onCluster.R**</a> and <a href="fitNullModel_onRAP.R">**fitNullModel_onRAP.R**</a> scripts designed for the HPC and RAP platforms, respectively. 

**HPC Cluster**
Use the <a href="HPC_slurm/generate_NullModel_args.sh">**HPC_slurm/generate_NullModel_args.sh**</a> script to generate the the corresponding arguments file `args_NullModel.txt`.
Then, submit the jobs by <a href="HPC_slurm/submit_NullModel.slurm">**HPC_slurm/submit_NullModel.slurm**</a> through the Slurm system using the following code in the terminal:
```
num_lines=$(wc -l < args_NullModel.txt)
sbatch --array=1-${num_lines} submit_NullModel.slurm
```

**RAP Platform**
Use the <a href="RAP_SwissArmyKnife/generate_NullModel_command.sh">**RAP_SwissArmyKnife/generate_NullModel_command.sh**</a> script to generate the corresponding execution commands file `SAK_batch_null_model.txt`.
Then, submit the job via Swiss Army Knife.


**Input:**
`phenotype_all.txt`. The output from Step1.1, including phenotype, covariates, and polygenic effects.
**Output:**
`objNull_chrA.rda`. Multiple RData files containing null model, where A represents the chromosome. Users can rename this file in <a href="RAP_SwissArmyKnife/generate_NullModel_command.sh">**RAP_SwissArmyKnife/generate_NullModel_command.sh**</a> or <a href="HPC_slurm/generate_NullModel_args.sh">**HPC_slurm/generate_NullModel_args.sh**</a>.


### Step 2a: Individual analysis for common varaints

Perform single-variant analysis for common and low-frequency variants across the genome using the SurvSTAAR pipeline.

<a href="HPC_script/IndividualAnalysisGDS_onCluster.R">**HPC_script/IndividualAnalysisGDS_onCluster.R**</a> and <a href="IndividualAnalysisGDS_onRAP.R">**IndividualAnalysisGDS_onRAP.R**</a> scripts designed for the aGDS format files to perform individual analysis.
<a href="HPC_script/IndividualAnalysisPlink_onCluster.R">**HPC_script/IndividualAnalysisPlink_onCluster.R**</a> and <a href="IndividualAnalysisPlink_onRAP.R">**IndividualAnalysisPlink_onRAP.R**</a> scripts designed for the plink format files to perform individual analysis.

First, use <a href="split_jobs/split_jobs_IndividualGDS.R">**split_jobs/split_jobs_IndividualGDS.R**</a> or <a href="split_jobs/split_jobs_IndividualPlink.R">**split_jobs/split_jobs_IndividualPlink.R**</a> to divide all variants into multiple sub-jobs.
**Output:** A numeric vector indicating the number of sub-jobs for each of the 22 chromosomes.

**HPC Cluster**
- Copy and paste the output numeric vector to set `chr_Nsets` in <a href="HPC_slurm/generate_IndividualGDS_args.sh">**HPC_slurm/generate_IndividualGDS_args.sh**</a> or <a href="HPC_slurm/generate_IndividualPlink_args.sh">**HPC_slurm/generate_IndividualPlink_args.sh**</a>, ensuring that the first element is 0.
- Run `sh generate_IndividualGDS_args.sh` or `sh generate_IndividualPlink_args.sh` in the terminal to generate the corresponding arguments file `args_IndividualGDS.txt` or `args_IndividualPlink.txt`.
- Submit the jobs by <a href="HPC_slurm/submit_IndividualGDS.slurm">**HPC_slurm/submit_IndividualGDS.slurm**</a> or <a href="HPC_slurm/submit_IndividualPlink.slurm">**HPC_slurm/submit_IndividualPlink.slurm**</a> through the Slurm system using the following code in the terminal:
  ```
  num_lines=$(wc -l < args_xxx.txt)
  sbatch --array=1-${num_lines} xxx.slurm
  ```
  ***NOTE***: The `args_xxx.txt` should be `args_IndividualGDS.txt` or `args_IndividualPlink.txt`, `xxx.slurm` should be `submit_IndividualGDS.slurm` or `submit_IndividualPlink.slurm`, depends on the file format.

**RAP Platform**
- Copy and paste the output numeric vector to set `chr_Nsets` in <a href="RAP_SwissArmyKnife/generate_IndividualGDS_command.sh">**RAP_SwissArmyKnife/generate_IndividualGDS_command.sh**</a> or <a href="RAP_SwissArmyKnife/generate_IndividualPlink_command.sh">**RAP_SwissArmyKnife/generate_IndividualPlink_command.sh**</a>, ensuring that the first element is 0.
- Run `sh generate_IndividualGDS_command.sh` or `sh generate_IndividualPlink_command.sh` in the terminal to generate the corresponding execution commands file `SAK_batch_individualGDS.txt` or `SAK_batch_individualPlink.txt`.
- Submit the job via Swiss Army Knife.  


**Output:**
`Individual_results_chrA_B.rda`. Multiple RData files containing individual analysis results, where A represents the chromosome and B represents the sub-job. Users can rename these files in <a href="HPC_slurm/generate_IndividualGDS_args.sh">**HPC_slurm/generate_IndividualGDS_args.sh**</a>, <a href="HPC_slurm/generate_IndividualPlink_args.sh">**HPC_slurm/generate_IndividualPlink_args.sh**</a>, <a href="RAP_SwissArmyKnife/generate_IndividualGDS_command.sh">**RAP_SwissArmyKnife/generate_IndividualGDS_command.sh**</a>, or <a href="RAP_SwissArmyKnife/generate_IndividualPlink_command.sh">**RAP_SwissArmyKnife/generate_IndividualPlink_command.sh**</a>.


### Step 2b: Gene-based test for rare variants

#### Step 2b.1 Gene-centric coding analysis
Perform gene-centric analysis for coding rare variants using the SurvSTAAR pipeline. The gene-centric coding analysis provides five functional categories to aggregate coding rare variants of each protein-coding gene: (1) putative loss of function (stop gain, stop loss, and splice) RVs, (2) missense RVs, (3) disruptive missense RVs, (4) putative loss of function and disruptive missense RVs, and (5) synonymous RVs.

<a href="HPC_script/GeneCentricCoding_onCluster.R">**HPC_script/GeneCentricCoding_onCluster.R**</a> and <a href="GeneCentricCoding_onRAP.R">**GeneCentricCoding_onRAP.R**</a> scripts designed for the aGDS format files to perform gene-centric coding analysis.

First, use <a href="split_jobs/split_jobs_coding_noncoding.R">**split_jobs/split_jobs_coding_noncoding.R**</a> to divide all genes into multiple sub-jobs.
**Output:** A numeric vector indicating the number of sub-jobs for each of the 22 chromosomes.


**HPC Cluster**
- Copy and paste the output numeric vector to set `chr_Nsets` in <a href="HPC_slurm/generate_Coding_args.sh">**HPC_slurm/generate_Coding_args.sh**</a>, ensuring that the first element is 0.
- Run `sh generate_Coding_args.sh` in the terminal to generate the corresponding arguments file `args_Coding.txt`.
- Submit the jobs by <a href="HPC_slurm/submit_Coding.slurm">**HPC_slurm/submit_Coding.slurm**</a> through the Slurm system using the following code in the terminal:
  ```
  num_lines=$(wc -l < args_Coding.txt)
  sbatch --array=1-${num_lines} submit_Coding.slurm
  ```

**RAP Platform**
- Copy and paste the output numeric vector to set `chr_Nsets` in <a href="RAP_SwissArmyKnife/generate_Coding_command.sh">**RAP_SwissArmyKnife/generate_Coding_command.sh**</a>, ensuring that the first element is 0.
- Run `sh generate_Coding_command.sh` in the terminal to generate the corresponding execution commands file `SAK_batch_coding.txt`.
- Submit the job via Swiss Army Knife.


**Output:**
`Coding_results_chrA_B.rda`. Multiple RData files containing gene-centric coding rare variants analysis results, where A represents the chromosome and B represents the sub-job. Users can rename these files in <a href="HPC_slurm/generate_Coding_args.sh">**HPC_slurm/generate_Coding_args.sh**</a> or <a href="RAP_SwissArmyKnife/generate_Coding_command.sh">**RAP_SwissArmyKnife/generate_Coding_command.sh**</a>.


#### Step 2b.2 Gene-centric noncoding analysis
Perform gene-centric analysis for noncoding rare variants using the the SurvSTAAR pipeline. The gene-centric noncoding analysis provides eight functional categories of regulatory regions to aggregate noncoding rare variants: (1) promoter RVs overlaid with CAGE sites, (2) promoter RVs overlaid with DHS sites, (3) enhancer RVs overlaid with CAGE sites, (4) enhancer RVs overlaid with DHS sites, (5) untranslated region (UTR) RVs, (6) upstream region RVs, (7) downstream region RVs, and (8) noncoding RNA (ncRNA) RVs.

<a href="HPC_script/GeneCentricNonCoding_onCluster.R">**HPC_script/GeneCentricNonCoding_onCluster.R**</a> and <a href="GeneCentricNonCoding_onRAP.R">**GeneCentricNonCoding_onRAP.R**</a> scripts designed for the aGDS format files to perform gene-centric noncoding analysis. 

<a href="HPC_script/ncRNA_onCluster.R">**HPC_script/ncRNA_onCluster.R**</a> and <a href="ncRNA_onRAP.R">**ncRNA_onRAP.R**</a> scripts designed for the aGDS format files to perform gene-centric noncoding analysis for ncRNA genes across the genome.


First, use <a href="split_jobs/split_jobs_coding_noncoding.R">**split_jobs/split_jobs_coding_noncoding.R**</a> and <a href="split_jobs/split_jobs_ncRNA.R">**split_jobs_ncRNA.R**</a> to divide all genes into multiple sub-jobs.
**Output:** A numeric vector indicating the number of sub-jobs for each of the 22 chromosomes.


**HPC Cluster**
- Copy and paste the output numeric vector separately to set `chr_Nsets` in <a href="HPC_slurm/generate_NonCoding_args.sh">**HPC_slurm/generate_NonCoding_args.sh**</a> and <a href="HPC_slurm/generate_ncRNA_args.sh">**HPC_slurm/generate_ncRNA_args.sh**</a>, ensuring that the first element is 0.
- Run `sh generate_NonCoding_args.sh` and `sh generate_ncRNA_args.sh` in the terminal to generate the corresponding arguments file `args_NonCoding.txt` and `args_ncRNA.txt`.
- Submit the jobs by <a href="HPC_slurm/submit_NonCoding.slurm">**HPC_slurm/submit_NonCoding.slurm**</a> and <a href="HPC_slurm/submit_ncRNA.slurm">**HPC_slurm/submit_ncRNA.slurm**</a> through the Slurm system using the following code in the terminal:
  ```
  num_lines=$(wc -l < args_NonCoding.txt)
  sbatch --array=1-${num_lines} submit_NonCoding.slurm
  
  num_lines=$(wc -l < args_ncRNA.txt)
  sbatch --array=1-${num_lines} submit_ncRNA.slurm
  ```

**RAP Platform**
- Copy and paste the output numeric vector separately to set `chr_Nsets` in <a href="RAP_SwissArmyKnife/generate_NonCoding_command.sh">**RAP_SwissArmyKnife/generate_NonCoding_command.sh**</a> and <a href="RAP_SwissArmyKnife/generate_ncRNA_command.sh">**RAP_SwissArmyKnife/generate_ncRNA_command.sh**</a>, ensuring that the first element is 0.
- Run `sh generate_NonCoding_command.sh` and `sh generate_ncRNA_command.sh` in the terminal to generate the corresponding execution commands file `SAK_batch_noncoding.txt` and `SAK_batch_ncRNA.txt`.
- Submit the job via Swiss Army Knife.


**Output:**
`NonCoding_results_chrA_B.rda` and `ncRNA_results_chrA_B.rda`. Multiple RData files containing gene-centric noncoding rare variants analysis results, where A represents the chromosome and B represents the sub-job. Users can rename these files in <a href="HPC_slurm/generate_NonCoding_args.sh">**HPC_slurm/generate_NonCoding_args.sh**</a>, <a href="RAP_SwissArmyKnife/generate_NonCoding_command.sh">**RAP_SwissArmyKnife/generate_NonCoding_command.sh**</a>, <a href="HPC_slurm/generate_ncRNA_args.sh">**HPC_slurm/generate_ncRNA_args.sh**</a>, or <a href="RAP_SwissArmyKnife/generate_ncRNA_command.sh">**RAP_SwissArmyKnife/generate_ncRNA_command.sh**</a>.


### Data processing and visualization

Use the provided script to organize the corresponding analysis results and generate Manhattan and QQ plots.
- <a href="SummaryFigures/Individual_AnalysisResults.R">**SummaryFigures/Individual_AnalysisResults.R**</a>
- <a href="SummaryFigures/Coding_AnalysisResults.R">**SummaryFigures/Coding_AnalysisResults.R**</a>
- <a href="SummaryFigures/NonCoding_AnalysisResults.R">**SummaryFigures/NonCoding_AnalysisResults.R**</a>
- <a href="SummaryFigures/ncRNA_AnalysisResults.R">**SummaryFigures/ncRNA_AnalysisResults.R**</a>

**Input:**
Revise the setwd() path to your results directory.
**Output:**
Each script generates an RData file containing all analysis results and two PNG files for the Manhattan and QQ plots.


## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
