# utility

## Comprehensive collection of Single-Cell Tumor-Infiltrating Lymphocyte Data

<img align="right" src="https://www.borch.dev/uploads/screpertoire/reference/figures/utility_hex.png" width="305" height="352">

### Introduction
The original intent of assembling a data set of publicly-available tumor-infiltrating T cells (TILs) with paired TCR sequencing was to expand 
and improve the [scRepertoire](https://github.com/ncborcherding/scRepertoire) R package. However, after some discussion, we decided to release 
the data set for everyone, a complete summary of the sequencing runs and the sample information can be found in the meta data of the Seurat object. 

An explanation of each variable is available [here](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/meta.data.headers.txt).

This involves several steps 1) loading the respective GE data, 2) harmonizing the data by sample and cohort information, 
3) iterating through automatic annotation, and 4) adding the TCR information. This information is stored in the meta data of the Seurat objects - an explanation of each variable is available [here](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/meta.data.headers.txt).

### Folder Structure
```
├── Data_conversion.Rmd
├── NEWS.txt
├── Processing_Utility.Rmd
├── README.md
├── Summarize_Data.Rmd
├── annotation
├── data
│   ├── SequencingRuns - 10x Outputs
│	  └── processedData - Processed .rds and larger combined cohorts
├── qc
└── summaryInfo
    ├── TcellSummaryTable.csv
    ├── cohortSummaryTable.csv
    ├── meta.data.headers.txt - what the meta data headers mean
    ├── sample.directory.xlsx - all the available data for the cohort
    ├── sessionInfo.txt - what I am running in terms of the pipeline
    └── tumorSummaryTable.csv
```

### Sample ID:

<img align="center" src="https://github.com/ncborcherding/utility/blob/dev/www/utility_info.png">

#### Cohort Information
Here is the current list of data sources, the number of cells that passed filtering by tissue type. Please cite the data if you are using utility!

|             | Blood | Juxta | LN   | Met | Normal | Tumor | Cancer Type | Date Added | Citation |
|-------------|-------|-------|------|-----|---|-------|-------------|------------|----------|
| CCR-20-4394 | 0     | 0     | 0    | 0   |0      | 26760 | Ovarian     | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33963000/) |
| EGAS00001004809| 0     | 0     | 0    | 0   | 0      | 181667 | Breast      | 3/30/22 |[cite](https://pubmed.ncbi.nlm.nih.gov/33958794/) |
| GSE114724   | 0     | 0     | 0    | 0   | 0      | 27651 | Breast      | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/29961579/) |
| GSE121636   | 12319 | 0     | 0    | 0   | 0      | 11436 | Renal       | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33504936/) |
| GSE123814   | 0     | 0     | 0    | 0   |0      | 77496 | Multiple    | 7/4/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/31359002/) |
| GSE139555   | 20664 | 0     | 0    | 0   | 69827  | 83301 | Multiple    | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/32103181/) |
| GSE145370   | 0     | 0     | 0    | 0   | 40916  | 66592 | Esophageal  | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33293583/) |
| GSE148190   | 6201  | 0     | 15644| 0   | 0      | 2263  | Melanoma    | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/32539073/) |
| GSE154826   | 0     | 0     | 0    | 0   | 13414   | 14491  | Lung    | 9/21/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/34767762/) |
| GSE159251   | 47721 | 0     | 5705 | 0   | 0      | 8355  | Melanoma    | 9/21/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/32539073/) |
| GSE162500   | 23401 | 3761  | 0    | 0   | 0      | 14644 | Lung        | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33514641/) |
| GSE164522   | 46027 | 0     | 46376|36648 | 86811 | 36990 | Colorectal | 6/25/22 | [cite](https://pubmed.ncbi.nlm.nih.gov/35303421/) |
| GSE176021   | 132673| 0     | 71062|32011 |128387 | 436608 | Lung      | 8/1/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/34290408/) |
| GSE179994   | 0     | 0     | 0    | 0   |0       | 140915 | Lung      | 3/30/22 |[cite](https://pubmed.ncbi.nlm.nih.gov/35121991/) |
| GSE180268   | 0     | 0     | 29699| 0   | 0      | 23215 | HNSCC      | 9/21/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/34471285/) |
| GSE181061   | 40429 | 0     | 0    | 0   | 27622  | 40429 | Renal      | 3/30/31 |[cite](https://pubmed.ncbi.nlm.nih.gov/35668194/) |
| GSE195486   | 0     | 0     | 0    | 0   | 0      | 122511 | Ovarian   | 6/25/22 |[cite](https://pubmed.ncbi.nlm.nih.gov/35427494/) |
| GSE200996   | 1211659| 0    | 0    | 0   | 0      | 86235  | HNSCC     | 7/15/22 | [cite](https://pubmed.ncbi.nlm.nih.gov/35803260/) | 
| GSE201425   | 27781 | 0    | 11350 |12253| 0      | 22888  | Biliary   | 1/18/23 | [cite](https://pubmed.ncbi.nlm.nih.gov/35982235/) | 
| GSE213243   | 18363 | 0    | 0     |0    | 0      | 22888  | Ovarian  | 1/18/23 | [cite](https://pubmed.ncbi.nlm.nih.gov/36248860/) | 
| PRJNA705465 | 30340 | 0     | 3505 | 0   | 15113  | 97966 | Renal      | 9/21/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33861994/) |

*****
### Methods

#### Single-Cell Data Processing
The filtered gene matrices output from Cell Ranger align function  from individual sequencing runs (10x Genomics, Pleasanton, CA) loaded into the R global environment. For each sequencing run cell barcodes were appended to contain a unique prefix to prevent issues with duplicate barcodes. The results were then ported into individual Seurat objects ([citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)), where the cells with > 10% mitochondrial genes and/or 2x natural log distribution of counts were excluded for quality control purposes. At the individual sequencing run level, doublets were estimated using the scDblFinder (v1.4.0) R package. All the sequencing runs across experiments were merged into a single Seurat Object using the merge() function. All the data was then normalized using the default settings and 2,000 variable genes were identified using the "vst" method. Next the data was scaled with the default settings and principal components were calculated for 40 components. 

#### Annotation of Cells

Automatic annotation was performed using the singler (v1.4.1) R package ([citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)) with the HPCA ([citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)) and DICE ([citation](https://pubmed.ncbi.nlm.nih.gov/30449622/)) data sets as references and the fine label discriminators. Individual sequencing runs were subsetted to run through the singleR algorithm in order to reduce memory demands. The output of all the singleR analyses were collated and appended to the meta data of the seurat object. Likewise, the ProjecTILs (v0.4.1) R Package ([citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)) was used for automatic annotation as a partially orthogonal approach. Consensus annotation was derived from all 3 databases (HPCA, Monaco, ProjecTILs) using a majority approach. No annotation designation was assigned to cells that returned NA for both singleR and ProjecTILs. 

#### Addition of TCR data

The filtered contig annotation T cell receptor (TCR) data for available sequencing runs were loaded into the R global environment. Individual contigs were combined using the combineTCR() function of scRepertoire (v2.0.0) R Package ([citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)). Clonotypes were assigned to barcodes and were multiple duplicate chains for individual cells were filtered to select for the top expressing contig by read count. The clonotype data was then added to the Seurat Object with proportion across individual patients being used to calculate frequency.

#### Session Info

Session Info for the initial data processing and analysis can be found [here](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/sessionInfo.txt).

*****
### Getting Data

Due to the size of the files, the  processed data outputs and code are available [here](https://zenodo.org/record/6325603).

*****
### Citations

As of right now, there is no citation associated with the assembled data set. However if using the data, please find the corresponding manuscript for 
each data set summarized above or can be found in the [summary table](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/cohortSummaryTable.csv). In addition, if using the processed data, feel free to modify the language in the 
methods section (above) and please cite the appropriate manuscripts of the software or references that were used.

#### Itemized List of the Software Used
* Seurat v4.1.0 - [citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)  
* singler v1.4.1 - [citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)  
* ProjecTILs v2.0.3 - [citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)
* UCell v1.0.0 - [citation](https://www.sciencedirect.com/science/article/pii/S2001037021002816?via%3Dihub)  
* scRepertoire v1.5.3 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)  

#### Itemized List of Reference Data Used
* Human Primary Cell Atlas (HPCA) - [citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)  
* Monaco Data Set (DICE) - [citation](https://pubmed.ncbi.nlm.nih.gov/30726743/)  

*****
### Future Directions

* Unified Dimensional Reduction of T cells with Cluster Annotations
* Data Hosting for Interactive Analysis
* Easy Submission Portal for Researchers to Add Data
* Using the Data to Build a Reference Atlas

There are areas in which we are actively hoping to develop to further facilitate the usefulness of the data set - if you have other suggestions, please reach out using the contact information below.

*****
### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

