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
* projecTILs v0.4.1 - [citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)
* UCell v1.0.0 - [citation](https://www.biorxiv.org/content/10.1101/2021.04.13.439670v1)  
* scRepertoire v1.3.2 - [citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)  

#### Itemized List of Reference Data Used
* Human Primary Cell Atlas (HPCA) - [citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)  
* Database Immune Cell Expression (DICE) - [citation](https://pubmed.ncbi.nlm.nih.gov/30449622/)  
* Immune-related Gene Sets - [citation](https://pubmed.ncbi.nlm.nih.gov/29961579/)

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 
