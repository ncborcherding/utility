# utility
Collection of Tumor-Infiltrating Lymphocyte Single-Cell Experiments with TCR sequencing data

### Introduction
The original intent of assembling a data set of publicly-available tumor-infiltrating T cells (TILs) with paired TCR sequencing was to expand 
and improve the [scRepertoire](https://github.com/ncborcherding/scRepertoire) R package. However, after some discussion, we decided to release 
the data set for everyone, a complete summary of the sequencing runs and the sample information can be found in the meta data of the Seurat object. 
This repository contains the fourth version of the data - with updates in adding data and the pipeline itself (v0.0.4). The major change here is that there is not currently a single integrated object - each sequencing run has a corresponding .rds and .h5ad. The meta data variables of each Seurat object is available [here](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/meta.data.headers.txt).

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
│	└── processedData - Processed .rds/.h5
├── qc
├── scGateDB
└── summaryInfo
    ├── TcellSummaryTable.csv
    ├── cohortSummaryTable.csv
    ├── meta.data.headers.txt - what the meta data headers mean
    ├── sample.directory.xlsx - all the available data for the cohort
    ├── sessionInfo.txt - what I am running in terms of the pipeline
    └── tumorSummaryTable.csv
```

#### Cohort Information
Here is the current list of data sources, the number of cells that passed filtering by tissue type. Please cite the data if you are using utility!

|             | Blood | Juxta | LN   | Met | Normal | Tumor | Cancer Type | Date Added | Citation |
|-------------|-------|-------|------|-----|---|-------|-------------|------------|----------|
| CCR-20-4394 | 0     | 0     | 0    | 0   |0      | 26760 | Ovarian     | 6/19/21 |[cite](https://clincancerres.aacrjournals.org/content/early/2021/06/10/1078-0432.CCR-20-4394) |
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
| GSE176021   | 132673| 0     | 71062|32011 |128387 | 436608 | Lung      | 8/1/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/34290408/) |
| GSE179994   | 0     | 0     | 0    | 0   |0       | 140915 | Lung      | 3/30/22 |[cite](https://pubmed.ncbi.nlm.nih.gov/35121991/) |
| GSE180268   | 0     | 0     | 29699| 0   | 0      | 23215 | HNSCC      | 9/21/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/34471285/) |
| GSE180268   | 40429 | 0     | 0    | 0   | 27622  | 40429 | Renal      | 3/30/31 |NA |
| PRJNA705465 | 30340 | 0     | 3505 | 0   | 15113  | 97966 | Renal      | 9/21/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33861994/) |

*****
### Methods

#### Single-Cell Data Processing
The filtered gene matrices output from Cell Ranger align function  from individual sequencing runs (10x Genomics, Pleasanton, CA) loaded into the R global environment. For each sequencing run cell barcodes were appended to contain a unique prefix to prevent issues with duplicate barcodes. The results were then ported into individual Seurat objects ([citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)), where the cells with > 10% mitochondrial genes and/or 2.5x natural log distribution of counts were excluded for quality control purposes. At the individual sequencing run level, doublets were estimated using the scDblFinder (v1.4.0) R package.  

#### Annotation of Cells

Automatic annotation was performed using the singler (v1.4.1) R package ([citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)) with the HPCA ([citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)) and Monaco ([citation](https://pubmed.ncbi.nlm.nih.gov/30726743/)) data sets as references and the fine label discriminators. Individual sequencing runs were subsetted to run through the singleR algorithm in order to reduce memory demands. The output of all the singleR analyses were collated and appended to the meta data of the seurat object. Likewise, the ProjecTILs (v2.0.3) R Package ([citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)) was used for automatic annotation as a partially orthogonal approach. Consensus annotation was derived from all 3 databases (HPCA, Monaco, ProjecTILs) using a majority approach. No annotation designation was assigned to cells that returned NA for both singleR and ProjecTILs. 

#### Addition of TCR data

The filtered contig annotation T cell receptor (TCR) data for available sequencing runs were loaded into the R global environment. Individual contigs were combined using the combineTCR() function of scRepertoire (v1.3.5) R Package ([citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)). Clonotypes were assigned to barcodes and were multiple duplicate chains for individual cells were filtered to select for the top expressing contig by read count. The clonotype data was then added to the Seurat Object with proportion across individual patients being used to calculate frequency.

#### Session Info

Session Info for the intial data processing and analysis can be found [here](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/sessionInfo.txt).

*****
### Getting Data

Due to the size of the files, the  processed data outputs and code are available [here](https://zenodo.org/record/6325603).

#### Most up-to-date version

#### Not a Static Resource 

Although the initial cohort of data is deposited and accessible - new data and analyses are being added. Check out the [dev branch](https://github.com/ncborcherding/utility/tree/dev) for changes and I am working through a number of other sequencing experiments, in order to provide the most up-to-date version (as I work), the following [link](https://drive.google.com/drive/folders/1Y8fGXIRxIfEk1BiQ4X2MC0CTznkXf_AW?usp=sharing) can be used.

*****
### Citations

As of right now, there is no citation associated with the assembled data set. However if using the data, please find the corresponding manuscript for 
each data set summarized above or can be found in the [summary table](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/cohortSummaryTable.csv). In addition, if using the processed data, feel free to modify the language in the 
methods section (above) and please cite the appropriate manuscripts of the software or references that were used.

#### Itemized List of the Software Used
* Seurat v4.1.p - [citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)  
* singler v1.4.1 - [citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)  
* ProjecTILs v2.0.3 - [citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)
* UCell v1.0.0 - [citation](https://www.sciencedirect.com/science/article/pii/S2001037021002816?via%3Dihub)  
* scRepertoire v1.3.5 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)  

#### Itemized List of Reference Data Used
* Human Primary Cell Atlas (HPCA) - [citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)  
* Monaco Data Set (Monaco) - [citation](https://pubmed.ncbi.nlm.nih.gov/30726743/)  

*****
### Future Directions

* Data Hosting for Interactive Analysis
* Easy Submission Portal for Researchers to Add Data
* Using the Data to Build a Reference Atlas

There are areas in which we are actively hoping to develop to further facilitate the usefulness of the data set - if you have other suggestions, please reach out using the contact information below.

*****
### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

******
### How to Support utility
This project has evolved into a bit of an undertaking, beyond time and effort, it has required caffeine. If you like to help, feel free to reach out or provide the caffeine directly:

[<img width="208" alt="snapshot-bmc-button" src="https://user-images.githubusercontent.com/22754118/175768906-2061f309-d038-4a0a-8f5b-5f250bf451cf.png">](https://www.buymeacoffee.com/theHumanBorch)
