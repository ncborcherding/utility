# utility
Collection of Tumor-Infiltrating Lymphocyte Single-Cell Experiments with TCR

### Introduction
The original intent of assembling a data set of publicly-available tumor-infiltrating T cells (TILs) with paired TCR sequencing was to expand 
and improve the [scRepertoire](https://github.com/ncborcherding/scRepertoire) R package. However, after some discussion, we decided to release 
the data set for everyone, a complete summary of the sequencing runs and the sample information can be found in the meta data of the Seurat object. 
This repository contains the code for the initial processing and annotating of the data set (we are calling this version 0.0.1). 
This involves several steps 1) loading the respective GE data, 2) harmonizing the data by sample and cohort information, 
3) iterating through automatic annotation, 4) unifying annotation via manual inspection and enrichment analysis, and 5) adding the TCR information. 

### Methods

#### Single-Cell Data Processing
The filtered gene matrices output from Cell Ranger align function  from individual sequencing runs  (10x Genomics, Pleasanton, CA) loaded into R global environment. For each sequencing run cell barcodes were appended to contain a unique prefix to prevent issues with duplicate barcodes. The reuslts were then ported into individual Seurat objects ([citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)), where the cells with > 10% mitochondrial genes and/or 2.5x natual log distribution of counts were exlcuded for quality control purposes. At the individual sequencing run level, doublets were estimated using the scDblFinder (v1.4.0) R package. All the sequencing runs across experiments were merged into a single Seurat Object using the merge() function. All the data was then normilized using the default settings and 2,000 variable genes were identified using the "vst" method. Next the data was scaled with the defualt settings and principal components were calculated for 40 components. Data was interegrated using the harmony (v1.0.0) R package ([citation](https://pubmed.ncbi.nlm.nih.gov/31740819/)) using both cohort and sample information to correct for batch effect with up to 20 iterations. The UMAP was created using the runUMAP() function in Seurat, using 20 dimensions of the harmony calculations. 

#### Annotation of Cells

Automatic annotation was performed using the singler (v1.4.1) R package (citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)) with the HPCA ([citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)) and DICE ([citation](https://pubmed.ncbi.nlm.nih.gov/30449622/)) data sets as references and the fine label discriminators. Individual sequencing runs were subsetted to run through the singleR algorithm in order to reduce memory demands. The output of all the singleR analyses were collated and appended to the meta data of the seurat object. Likewise, the ProjecTILs (v0.4.1) R Package ([citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)) was used for autmatic annotation as a partially orthagonal approach. Consensus annotation was derived from all 3 databases (HPCA, DICE, ProjecTILs) using a majority approach. No annotation designation was assigned to cells that returned NA for both singleR and ProjecTILs. Mixed annotations were designated with SingleR identified non-Tcells and ProjecTILs identified T cells. Cell type designations with less than 100 cells in the entire cohort were reduced to "other".

#### Addition of TCR data

### Getting Data

Due to the size of the files, the processed 10x Genomics outputs and processed data outputs from the anlaytical pipeline are available here (link to come).

### Citations

As of right now, there is no citation associated with the assembled data set. However if using the data, please find the corresponding manuscipt for 
each data set in the meta.data of the the single-cell object. In addition, if using the processed data, feel free to modify the language in the 
methods section (above) and please cite the appropriate manuscripts of the software or references that were used.

#### Itemized List of the Software Used
* Seurat v4.0.3 - [citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)  
* harmony v1.0 - [citation](https://pubmed.ncbi.nlm.nih.gov/31740819/)  
* singler v1.4.1 - [citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)  
* ProjecTILs v0.4.1 - [citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)
* UCell v1.0.0 - [citation](https://www.biorxiv.org/content/10.1101/2021.04.13.439670v1)  
* scRepertoire v1.3.2 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)  

#### Itemized List of Reference Data Used
* Human Primary Cell Atlas (HPCA) - [citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)  
* Database Immune Cell Expression (DICE) - [citation](https://pubmed.ncbi.nlm.nih.gov/30449622/)  
* Immune-related Gene Sets - [citation](https://pubmed.ncbi.nlm.nih.gov/29961579/)

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 
