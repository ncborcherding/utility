# uTILity

## Comprehensive collection of Single-Cell Tumor-Infiltrating Lymphocyte Data

<img align="right" src="https://github.com/ncborcherding/utility/blob/main/www/utility_hex.png" width="305" height="352">

### Introduction
The original intent of assembling a data set of publicly-available tumor-infiltrating T cells (TILs) with paired TCR sequencing was to expand 
and improve the [scRepertoire](https://github.com/ncborcherding/scRepertoire) R package. However, after some discussion, we decided to release 
the data set for everyone, a complete summary of the sequencing runs and the sample information can be found in the meta data of the Seurat object. 

This involves several steps 1) loading the respective GE data, 2) harmonizing the data by sample and cohort information, 3) iterating through automatic annotation, 
and 4) adding the TCR information. This information is stored in the meta data of the Seurat objects - 
an explanation of each variable is available [here](https://github.com/ncborcherding/utility/blob/dev/summaryInfo/meta.data.headers.txt).


### Folder Structure
```
├── Data_conversion.Rmd
├── NEWS.txt
├── Processing_Utility.Rmd
├── README.md
├── Summarize_Data.Rmd
├── data
│   ├── SequencingRuns - 10x Outputs
│	└── processedData - Processed .rds and larger combined cohorts
├── outputs
│   └── qc - plots for quality control purposes
└── summaryInfo
    ├── TcellSummaryTable.csv
    ├── cohortSummaryTable.csv
    ├── meta.data.headers.txt - what the meta data headers mean
    ├── sample.directory.xlsx - all the available data for the cohort
    ├── sessionInfo.txt - what I am running in terms of the pipeline
    └── tumorSummaryTable.csv
```

### Sample ID:

<img align="center" src="https://github.com/ncborcherding/utility/blob/main/www/utility_info.png">


#### Cohort Information
Here is the current list of data sources, the number of cells that passed filtering by tissue type. **Please cite** the data if you are using uTILity.


|                  | Tumor  | Normal | Blood  | Juxta | LN    | Met   | Cancer Type   | Citations                                         |
|------------------|--------|--------|--------|-------|-------|-------|---------------|---------------------------------------------------|
| CCR-20-4394      | 26760  | 0      | 0      | 0     | 0     | 0     | Ovarian       | [cite](https://pubmed.ncbi.nlm.nih.gov/33963000/) |
| EGAS00001004809  | 181667 | 0      | 0      | 0     | 0     | 0     | Breast        | [cite](https://pubmed.ncbi.nlm.nih.gov/33958794/) |
| GSE114724        | 27651  | 0      | 0      | 0     | 0     | 0     | Breast        | [cite](https://pubmed.ncbi.nlm.nih.gov/29961579/) |
| GSE121636        | 11436  | 0      | 12319  | 0     | 0     | 0     | Renal         | [cite](https://pubmed.ncbi.nlm.nih.gov/33504936/) |
| GSE123814        | 78034  | 0      | 0      | 0     | 0     | 0     | Multiple      | [cite](https://pubmed.ncbi.nlm.nih.gov/31359002/) |
| GSE139555        | 93160  | 78625  | 25363  | 0     | 0     | 0     | Multiple      | [cite](https://pubmed.ncbi.nlm.nih.gov/32103181/) |
| GSE145370        | 66592  | 40916  | 0      | 0     | 0     | 0     | Esophageal    | [cite](https://pubmed.ncbi.nlm.nih.gov/33293583/) |
| GSE148190        | 2263   | 0      | 6201   | 0     | 15644 | 0     | Melanoma      | [cite](https://pubmed.ncbi.nlm.nih.gov/32539073/) |
| GSE154826        | 14491  | 13414  | 0      | 0     | 0     | 0     | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/34767762/) |
| GSE159251        | 8356   | 0      | 47721  | 0     | 5705  | 0     | Melanoma      | [cite](https://pubmed.ncbi.nlm.nih.gov/32539073/) |
| GSE162500        | 14644  | 0      | 23401  | 3761  | 0     | 0     | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/33514641/) |
| GSE164522        | 36990  | 86811  | 46027  | 0     | 46376 | 36648 | Colorectal    | [cite](https://pubmed.ncbi.nlm.nih.gov/35303421/) |
| GSE168844        | 0      | 0      | 55302  | 0     | 0     | 0     | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/36219677/) |
| GSE176021        | 436609 | 128411 | 132673 | 0     | 71063 | 32011 | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/34290408/) |
| GSE179994        | 78574  | 0      | 0      | 0     | 0     | 62341 | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/35121991/) |
| GSE180268        | 23215  | 0      | 0      | 0     | 29699 | 0     | HNSCC         | [cite](https://pubmed.ncbi.nlm.nih.gov/34471285/) |
| GSE181061        | 40429  | 27622  | 37426  | 0     | 0     | 0     | Renal         | [cite](https://pubmed.ncbi.nlm.nih.gov/35668194/) |
| GSE185206        | 163294 | 17231  | 0      | 0     | 9820  | 0     | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/37001526/) |
| GSE195486        | 122512 | 0      | 0      | 0     | 0     | 0     | Ovarian       | [cite](https://pubmed.ncbi.nlm.nih.gov/35427494/) |
| GSE200218        | 0     | 0       | 0      | 0     | 0     | 18495 | Melanoma      | [cite](https://pubmed.ncbi.nlm.nih.gov/35803246/) |
| GSE200996        | 86235 | 0       | 152722 | 0     | 0     | 0     | HNSCC         | [cite](https://pubmed.ncbi.nlm.nih.gov/35803260/) |
| GSE201425        | 22888 | 0       | 27781  | 0     | 11350 | 12253 | Biliary       | [cite](https://pubmed.ncbi.nlm.nih.gov/35982235/) | 
| GSE211504        | 0     | 0       | 33685  | 0     | 0     | 0     | Melanoma      | [cite](https://pubmed.ncbi.nlm.nih.gov/35907015/) |
| GSE212217        | 0     | 0       | 229505 | 0     | 0     | 0     | Endometrial   | [cite](https://pubmed.ncbi.nlm.nih.gov/36301137/) |
| GSE213243        | 2835  | 0       | 18363  | 0     | 0     | 2693  | Ovarian       | [cite](https://pubmed.ncbi.nlm.nih.gov/36248860/) |
| GSE215219        | 26303 | 0       | 66000  | 0     | 0     | 0     | Lung          | [cite](https://pubmed.ncbi.nlm.nih.gov/37476074/) |
| GSE227708        | 53087 | 0       | 0      | 0     | 0     | 0     | Merkel Cell   | [cite](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) |
| GSE242477        | 41595 | 0       | 21595  | 0     | 0     | 0     | Melanoma      | [cite](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242477) |
| PRJNA705464      | 98892 | 15113   | 30340  | 0     | 3505  | 0     | Renal         | [cite](https://pubmed.ncbi.nlm.nih.gov/33861994/) |

*****
### Methods

#### Single-Cell Data Processing
The filtered gene matrices output from Cell Ranger align function  from individual sequencing runs (10x Genomics, Pleasanton, CA) loaded into the R global environment. For each sequencing run cell barcodes were appended to contain a unique prefix to prevent issues with duplicate barcodes. The results were then ported into individual Seurat objects ([citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)), where the cells with > 10% mitochondrial genes and/or 2.5x standard deviation from the mean of features were excluded for quality control purposes. At the individual sequencing run level, doublets were estimated using the scDblFinder (v1.4.0) R package. 

#### Annotation of Cells

Automatic annotation was performed using the singler (v2.2.0) R package ([citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)) with the HPCA ([citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)) and DICE ([citation](https://pubmed.ncbi.nlm.nih.gov/30449622/)) data sets as references and the fine label discriminators. Individual sequencing runs were subsetted to run through the singleR algorithm in order to reduce memory demands. The output of all the singleR analyses were collated and appended to the meta data of the seurat object. Likewise, the Azimuth (v0.4.6.9004) R Package ([citation](https://pubmed.ncbi.nlm.nih.gov/34062119/) was used for automatic annotation as a partially orthogonal approach. 

#### Addition of TCR data

The filtered contig annotation T cell receptor (TCR) data for available sequencing runs were loaded into the R global environment. Individual contigs were combined using the combineTCR() function of scRepertoire (v2.0.0) R Package ([citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)). Clonotypes were assigned to barcodes and were multiple duplicate chains for individual cells were filtered to select for the top expressing contig by read count. The clonotype data was then added to the Seurat Object with proportion across individual patients being used to calculate frequency.

#### Session Info


Session Info for the initial data processing and analysis can be found [here](https://github.com/ncborcherding/utility/blob/main/summaryInfo/sessionInfo.txt).

*****
### Citations

As of right now, there is no citation associated with the assembled data set. However if using the data, please find the corresponding manuscript for 
each data set summarized above or can be found in the [summary table](https://github.com/ncborcherding/utility/blob/main/summaryInfo/cohortSummaryTable.csv). In addition, if using the processed data, feel free to modify the language in the methods section (above) and please cite the appropriate manuscripts of the software or references that were used.

#### Itemized List of the Software Used
* Seurat v5.0.0 - [citation](https://pubmed.ncbi.nlm.nih.gov/37231261/)  
* singler v2.2.0 - [citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)  
* Azimuth v0.4.6.9004 - [citation](https://pubmed.ncbi.nlm.nih.gov/34062119/) 
* scRepertoire v2.0.0 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)  

#### Itemized List of Reference Data Used
* Human Primary Cell Atlas (HPCA) - [citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)  
* Monaco Data Set (Monaco) - [citation](https://pubmed.ncbi.nlm.nih.gov/30726743/)
* PBMC reference - [citation](https://pubmed.ncbi.nlm.nih.gov/31178118/)

*****
### Future Directions

* Unified Dimensional Reduction of T cells with Cluster Annotations
* Data Hosting for Interactive Analysis
* Easy Submission Portal for Researchers to Add Data
* Using the Data to Build a Reference Atlas

There are areas in which we are actively hoping to develop to further facilitate the usefulness of the data set - if you have other suggestions, please reach out using the contact information below.

*****
### License

The data and analysis of uTILity is provided under a CC BY-ND 4.0 license, please feel free to remix, transform, and build upon the material. However, the intent of this resource is noncommercial, if using the data as a nonacademic institution, you are in violation of the lisence agreement. Please find out more information [here](https://github.com/ncborcherding/utility/blob/main/LICENSE.txt).


*****
### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

