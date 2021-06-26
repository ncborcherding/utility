# utility
Collection of Tumor-Infiltrating Lymphocyte Single-Cell Experiments with TCR sequencing data

### Introduction
The original intent of assembling a data set of publicly-available tumor-infiltrating T cells (TILs) with paired TCR sequencing was to expand 
and improve the [scRepertoire](https://github.com/ncborcherding/scRepertoire) R package. However, after some discussion, we decided to release 
the data set for everyone, a complete summary of the sequencing runs and the sample information can be found in the meta data of the Seurat object. 
This repository contains the code for the initial processing and annotating of the data set (we are calling this version 0.0.1). 
This involves several steps 1) loading the respective GE data, 2) harmonizing the data by sample and cohort information, 
3) iterating through automatic annotation, 4) unifying annotation via manual inspection and enrichment analysis, and 5) adding the TCR information. This information is stored in the meta data of the Seurat object - an explanation of each variable is available [here](https://github.com/ncborcherding/utility/blob/main/meta.data.headers.txt).

#### Cohort Information
Here is the current list of data sources, the number of cells that passed filtering by tissue type. Please cite the data if you are using utility!

|             | Blood | Juxta | LN   | Normal | Tumor | Cancer Type | Date Added | Citation |
|-------------|-------|-------|------|--------|-------|-------------|------------|----------|
| CCR-20-4394 | 0     | 0     | 0    | 0      | 26767 | Ovarian     | 6/19/21 |[cite](https://clincancerres.aacrjournals.org/content/early/2021/06/10/1078-0432.CCR-20-4394) |
| GSE114724   | 0     | 0     | 0    | 0      | 27651 | Breast      | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/29961579/) |
| GSE121636   | 12322 | 0     | 0    | 0      | 11436 | Renal       | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33504936/) |
| GSE139555   | 20676 | 0     | 0    | 70114  | 83336 | Multiple    | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/32103181/) |
| GSE144469   | 0     | 0     | 0    | 41311  | 86436 | Colorectal  | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/32603654/) |
| GSE145370   | 0     | 0     | 0    | 41081  | 66739 | Esophageal  | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33293583/) |
| GSE148190   | 0     | 0     | 6370 | 0      | 2263  | Melanoma    | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/32539073/) |
| GSE162500   | 23426 | 3771  | 0    | 0      | 14647 | Lung        | 6/19/21 |[cite](https://pubmed.ncbi.nlm.nih.gov/33514641/) |

*****
### Methods

#### Single-Cell Data Processing
The filtered gene matrices output from Cell Ranger align function  from individual sequencing runs (10x Genomics, Pleasanton, CA) loaded into the R global environment. For each sequencing run cell barcodes were appended to contain a unique prefix to prevent issues with duplicate barcodes. The results were then ported into individual Seurat objects ([citation](https://pubmed.ncbi.nlm.nih.gov/34062119/)), where the cells with > 10% mitochondrial genes and/or 2.5x natural log distribution of counts were excluded for quality control purposes. At the individual sequencing run level, doublets were estimated using the scDblFinder (v1.4.0) R package. All the sequencing runs across experiments were merged into a single Seurat Object using the merge() function. All the data was then normalized using the default settings and 2,000 variable genes were identified using the "vst" method. Next the data was scaled with the default settings and principal components were calculated for 40 components. Data was integrated using the harmony (v1.0.0) R package ([citation](https://pubmed.ncbi.nlm.nih.gov/31740819/)) using both cohort and sample information to correct for batch effect with up to 20 iterations. The UMAP was created using the runUMAP() function in Seurat, using 20 dimensions of the harmony calculations. 

#### Annotation of Cells

Automatic annotation was performed using the singler (v1.4.1) R package ([citation](https://pubmed.ncbi.nlm.nih.gov/30643263/)) with the HPCA ([citation](https://pubmed.ncbi.nlm.nih.gov/24053356/)) and DICE ([citation](https://pubmed.ncbi.nlm.nih.gov/30449622/)) data sets as references and the fine label discriminators. Individual sequencing runs were subsetted to run through the singleR algorithm in order to reduce memory demands. The output of all the singleR analyses were collated and appended to the meta data of the seurat object. Likewise, the ProjecTILs (v0.4.1) R Package ([citation](https://pubmed.ncbi.nlm.nih.gov/34017005/)) was used for automatic annotation as a partially orthogonal approach. Consensus annotation was derived from all 3 databases (HPCA, DICE, ProjecTILs) using a majority approach. No annotation designation was assigned to cells that returned NA for both singleR and ProjecTILs. Mixed annotations were designated with SingleR identified non-Tcells and ProjecTILs identified T cells. Cell type designations with less than 100 cells in the entire cohort were reduced to "other". Automated annotations were checked manually using canonical marker genes and gene enrichment analysis performed using UCell (v1.0.0) R package ([citation](https://www.biorxiv.org/content/10.1101/2021.04.13.439670v1)).

#### Addition of TCR data

The filtered contig annotation T cell receptor (TCR) data for available sequencing runs were loaded into the R global environment. Individual contigs were combined using the combineTCR() function of scRepertoire (v1.3.2) R Package ([citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400693/)). Clonotypes were assigned to barcodes and were multiple duplicate chains for individual cells were filtered to select for the top expressing contig by read count. The clonotype data was then added to the Seurat Object with proportion across individual patients being used to calculate frequency.

#### Session Info

Session Info for the intial data processing and analysis can be found [here](https://github.com/ncborcherding/utility/blob/main/sessionInfo.txt).

*****
### Getting Data

Due to the size of the files, the processed 10x Genomics outputs and processed data outputs from the analytical pipeline are available [here](https://zenodo.org/record/4995299).

******
### Working with the Data

There is no assigned clusters across the single-cells. In order to incorporate additional data sets for use, we opted not to cluster. Clustering can easily be performed using the code below:

```
seuratObject <- FindNeighbors(seuratObject, reduction = "harmony") 
seuratObject <- FindClusters(seuratObject)
```

If subsetting the data, either by cell type or other variables in the meta data, re-intergration using Harmony will be required:

```
subset.1 <- subset(seuratObject, "Consensus.Major" == "T_cell")
subset.1 <- RunHarmony(subset.1, group.by.vars = c("Cohort", "Sample"), max.iter.harmony = 20)
subset.1 <- RunUMAP(subset.1, reduction = "harmony", dims = 1:20)
```

*****
### Citations

As of right now, there is no citation associated with the assembled data set. However if using the data, please find the corresponding manuscript for 
each data set in the meta.data of the single-cell object or the [summary table](https://github.com/ncborcherding/utility/blob/main/cohortSummaryTable.txt). In addition, if using the processed data, feel free to modify the language in the 
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

*****
### Future Directions

* Data Hosting for Interactive Analysis
* Easy Submission Portal for Researchers to Add Data
* Using the Data to Build a Reference Atlas

There are areas in which we are actively hoping to develop to further facilitate the usefulness of the data set - if you have other suggestions, please reach out using the contact information below.

*****
### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 
