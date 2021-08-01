#### Script for adding additional data sets to utility

#Some assumptions need to be made:
#1) the processed data outputs appear ./data/SequencingRuns
#2) Data information has been added to the sample directory
#3) The sample information provided is unique

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(scDblFinder))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(ProjecTILs))
suppressPackageStartupMessages(library(scRepertoire))
suppressPackageStartupMessages(library(UCell))

ref <- load.reference.map("./annotation/ref_TILAtlas_mouse_v1.rds")
HPCA <- HumanPrimaryCellAtlasData()
DICE <- DatabaseImmuneCellExpressionData()
signature.list <- readRDS("./data/processedData/signature.list.rds")

#Load most-up-to-date seurat object
meta <- readRDS(file = "./data/ProcessedData/meta.rds")

'%!in%' <- Negate('%in%')
file_list <- list.files("./data/SequencingRuns")
new.file_list <- file_list[file_list %!in% unique(meta$orig.ident)]
new.file_list <- new.file_list[!grepl("Icon", new.file_list)]

#for (y in seq_along(new.file_list)){
for (y in seq_along(new.file_list)){
  tmp <-  Read10X(paste0("./data/SequencingRuns/", new.file_list[y]))
  
  ##Several data sets do not have the MT in front of the mitochondria genes
  mito <- c("ATP8", "ATP6", "CO1", "CO2", "CO3", "CYB", "ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6", "RNR2", "TA", "TR", "TN", "TD", "TC", "TE", "TQ", "TG", "TH", "TI", "TL1", "TL2", "TK", "TM", "TF", "TP", "TS1", "TS2", "TT", "TW", "TY", "TV", "RNR1")
  
  x <- which(rownames(tmp) %in% mito)
  if (length(x) > 0) {
    z <- rownames(tmp)[x] 
    z<- paste0("MT-", z)
    rownames(tmp)[x] <- z
  }
  tmp <- CreateSeuratObject(counts = tmp)
  tmp <- subset(tmp, subset = nFeature_RNA > 100) #added filter to remove 0s and too sparse cells
  tmp  <- RenameCells(object = tmp , new.names = paste0(new.file_list[y], "_", rownames(tmp[[]])))
  
  tmp[["mito.genes"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  
  p1 <- VlnPlot(object = tmp, features = c("nCount_RNA")) + theme(legend.position = "none")
  p2 <- VlnPlot(object = tmp, features = c("nFeature_RNA")) + theme(legend.position = "none")
  p3 <- VlnPlot(object = tmp, features = c("mito.genes")) + theme(legend.position = "none")
  
  pdf(paste0("./qc/", new.file_list[y], ".pdf"), height = 8, width=12)
  grid.arrange(p1, p2, p3, ncol = 3)
  dev.off()
  
  ###########################
  #Here is the filtering step
  ############################
  standev <- sd(log(tmp$nFeature_RNA))*2.5 #cutting off above standard deviation of 2.5
  mean <- mean(log(tmp$nFeature_RNA))
  cut <- round(exp(standev+mean))
  tmp <- subset(tmp, subset = mito.genes < 10 & nFeature_RNA < cut)
  
  ###########################################
  #Estimate Doublets for Each Sequencing Run
  ############################################
  sce <- as.SingleCellExperiment(tmp)
  sce <- scDblFinder(sce, BPPARAM=MulticoreParam(3))
  doublets <- data.frame(db.weight.score = sce$scDblFinder.score, db.ratio = sce$scDblFinder.weighted, 
                         db.class = sce$scDblFinder.class, db.score = sce$scDblFinder.score)
  rownames(doublets) <- rownames(sce@colData)
  tmp <- AddMetaData(tmp, doublets)
  rm(sce)
  
  ####Adding meta data
  directory <- readxl::read_xlsx("./summaryInfo/sample.directory.xlsx") #Meta.data
  label <- stringr::str_split(rownames(tmp[[]]), "_", simplify = T)[,1]
  tmp[["orig.ident"]] <- label
  meta <- tmp[[]]
  rownames <- rownames(meta)
  meta <- merge(meta, directory, by.x = "orig.ident", by.y = "SampleLabel")
  meta <- meta[,9:15]
  rownames(meta) <- rownames
  
  tmp <- AddMetaData(tmp, meta)
  
  
  ####################
  #Adding Annotation
  #######################
  
  #Singler
  tmp.2 <- tmp@assays[["RNA"]]@counts
  ####This approach for matrix conversion saves some memory
  tmp.2 <- tmp.2[tabulate(summary(tmp.2)$i) != 0, , drop = FALSE]
  tmp.2 <- as.matrix(tmp.2)
  com.res1 <- SingleR(tmp.2, ref=HPCA, labels=HPCA$label.fine, assay.type.test=1)
  saveRDS(com.res1, file = paste0("./annotation/singler/",  new.file_list[y], "_HPCA.singler_output.rds"))
  com.res2 <- SingleR(tmp.2, ref=DICE, labels=DICE$label.fine, assay.type.test=1)
  saveRDS(com.res2, file = paste0("./annotation/singler/",  new.file_list[y], "_DICE.singler_output.rds"))
  rm(tmp.2)

  df <- data.frame("HPCA.first.labels" = com.res1$first.labels, "HPCA.labels" = com.res1$labels, "HPCA.pruned.labels" = com.res1$pruned.labels, 
                       "DICE.first.labels" = com.res2$first.labels, "DICE.labels" = com.res2$labels, "DICE.pruned.labels" = com.res2$pruned.labels)
  rownames(df) <- rownames(com.res1)
  tmp <- AddMetaData(tmp,  df)
  
  #Projectil
  query.projected <- make.projection(tmp, ref = ref, ncores = 1)
  query.projected <- cellstate.predict(ref = ref, query = query.projected)
  meta <- query.projected[[c("functional.cluster", "functional.cluster.conf")]]
  rownames(meta ) <- stringr::str_remove(rownames(meta), "Q_")
  tmp <- AddMetaData(tmp, meta)
  saveRDS(meta, file = paste0("./annotation/projectil/",  new.file_list[y],  "_output.rds"))
  rm(query.projected)
  
  #####
  #consensus Annotation
  ####
  consensus.df <- data.frame(DICE = stringr::str_split(tmp[[]]$DICE.pruned.labels, ",", simplify = TRUE)[,1], 
                             HPCA = stringr::str_split(tmp[[]]$HPCA.pruned.labels, ":", simplify = TRUE)[,1], 
                             PTIL = tmp[[]]$functional.cluster)
  
  
  consensus.df$DICE<- gsub(" ", "_", consensus.df$DICE)
  consensus.df$HPCA<- paste0(consensus.df$HPCA, "s")
  
  
  consensus.df$PTIL <- ifelse(!is.na(consensus.df$PTIL), "T_cells", NA)
  consensus.major <- NULL
  for (i in 1:nrow(consensus.df)) {
    T.length <- length(which(consensus.df[i,] %in% "T_cells"))
    
    if (T.length >= 2) {
      consensus.major[i] <- "T_cell"
    } else {
      cellType <- unlist(consensus.df[i,])
      if (length(which(is.na(cellType))) == 3) {
        consensus.major[i] <- "no.annotation"
        next()
      }
      if (cellType[1] == cellType[2] & is.na(cellType[3]) | length(which(is.na(cellType))) == 2) {
        consensus.major[i] <- cellType[1]
      } else {
        consensus.major[i] <- "mixed.annotation"
      }
    }
  }
  consensus.major <- data.frame(consensus.major)
  rownames(consensus.major) <- rownames(tmp[[]])
  
  consensus.df <- data.frame(DICE = tmp[[]]$DICE.labels,  
                             HPCA = tmp[[]]$HPCA.pruned.labels,
                             PTIL = tmp[[]]$functional.cluster)
  
  
  consensus.Tcell <- NULL
  for (i in 1:nrow(consensus.df)) {
    CD8.length <- length(which(grepl("CD8", consensus.df[i,])))
    CD4.length <- length(which(grepl("CD4", consensus.df[i,])))
    Treg.length <- length(which(grepl("Treg|TREG", consensus.df[i,])))
    if (CD8.length >= 2){
      consensus.Tcell[i] <- "CD8"
    } else if (CD4.length >= 2) {
      consensus.Tcell[i] <- "CD4"
    } else if (Treg.length >= 2) {
      consensus.Tcell[i] <- "Treg"
    } else {
      consensus.Tcell[i] <- NA
    }
  }
  consensus.Tcell <- data.frame(consensus.Tcell)
  rownames(consensus.Tcell) <- rownames(tmp[[]])
  
  tmp <- AddMetaData(tmp, consensus.major)
  tmp <- AddMetaData(tmp, consensus.Tcell)
  

  
  scores <- ScoreSignatures_UCell(tmp@assays$RNA@data, features=signature.list)
  tmp <- AddMetaData(tmp, as.data.frame(scores))
    
  if (y == 1) {
    list <- tmp
  } else {
    list <- merge(x=list, y=tmp)
  }
  rm(tmp)
}

saveRDS(list, file = "tmp.rds")
rm(ref)
rm(HPCA)
rm(DICE)
rm(meta)


list <- NormalizeData(list, verbose = FALSE, assay = "RNA")
list <- FindVariableFeatures(list, selection.method = "vst", 
                             nfeatures = 2000, verbose = FALSE, assay = "RNA")
list <- ScaleData(object = list, verbose = FALSE)
list <- RunPCA(object = list, npcs = 40, verbose = FALSE)

#######################
#Adding the meta data 
#######################


cohortSummary <- table(list$Cohort, list$Type)
write.csv(cohortSummary, file = "./summaryInfo/cohortSummaryTable.csv")

####################################################
#Correcting for Cohort and Sample and getting UMAP
###################################################
library(harmony)
list <- RunHarmony(list, group.by.vars = c("Cohort", "Sample"), max.iter.harmony = 20)


list <- RunUMAP(list, reduction = "harmony", dims = 1:20)
saveRDS(list, file = "./data/ProcessedData/filtered_seuratObjects_harmony.rds")

#Simplifying annotation
table <- table(list[["consensus.major"]])
table <- table[table > 100]
list$consensus.major <- ifelse(list$consensus.major %in% names(table), list$consensus.major, "other")
list$consensus.major <- ifelse(list$consensus.major == "no.annotation", "other", list$consensus.major)
saveRDS(list, file = "./data/ProcessedData/filtered_seuratObjects_harmony.rds")

mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(list$Tissue)))

DimPlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 20000)], group.by = "Tissue") + 
  scale_color_manual(values = mycolors) + 
  theme(plot.title = element_blank())
ggsave("./UMAP/TumorType.pdf", height = 3.5, width = 4.5)


mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(list$Type)))

DimPlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 20000)], group.by = "Type") + 
  scale_color_manual(values = mycolors) + 
  theme(plot.title = element_blank())
ggsave("./UMAP/TissueType.pdf", height = 3.5, width = 4.25)

mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(15)
DimPlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 20000)], group.by = "db.class") + 
  scale_color_manual(values = mycolors[c(1,3)]) + 
  theme(plot.title = element_blank())
ggsave("./UMAP/scDoublet.pdf", height = 3.5, width = 4.5)


#Examining the annotations using the sample function to reduce load - graphing 20,000 cell instead of 500,000+

library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(list$consensus.major)))

set.seed(123)
x <- sample(nrow(list[[]]), 20000)

dir.create("./UMAP")

DimPlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 20000)], group.by = "consensus.major") + 
  scale_color_manual(values = mycolors) + 
  theme(plot.title = element_blank())
ggsave("./UMAP/major.consensus.pdf", height = 3.5, width = 5.5)

DimPlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 20000)], group.by = "consensus.Tcell") + 
  scale_color_manual(values = mycolors[c(1,3)], na.value="grey") + 
  theme(plot.title = element_blank())
ggsave("./UMAP/Tcell.consensus.pdf", height = 3.5, width = 4)


# Looking at canonical marker expression

dir.create("DataAnalysis/UMAP/LineageMarkers")

file_list <- list.files("./data/ProcessedData/markers.genes")
file_list <- file_list[grepl(".txt", file_list)]
files <- file.path(paste0("./data/processedData/markers.genes/", file_list))

marker_list <- list()
for (i in 1:length(files)) {
  marker_list[[i]] <- read.delim(paste0(files[i]), col.names = FALSE)
  for (j in seq_len(nrow(marker_list[[i]]))) {
    marker_list[[i]][j,] <- toupper( marker_list[[i]][j,])
  }
}
names <- stringr::str_remove(file_list, ".txt")
names(marker_list) <- names

################################
#Graphing canonical marker genes
#################################
#I like to use the schex package in order to look at the percentage of cells in the area of the UMAP with GeneX expression
#This prevents the issue with overlap that can be seen in DotPlots, reduces the size of the subsequent visualization files, and as UMAP is based on nearest neighbor it make more sense to incorporate neighbors while visualizing (at least to me, I have a reviewer or two that have not been so inclined)

suppressPackageStartupMessages(library(schex))

list <- make_hexbin(list, 128, dimension_reduction = "UMAP")

DefaultAssay(list) <- "RNA" #This is probably not necessary as there is only one assay but it is an important step if the seurat object has integrated data assay


dir.create("./UMAP/LineageMarkers")
for (i in seq_along(marker_list)) {
  tmp <- as.character(unlist(marker_list[i]))
  for (j in seq_along(tmp)) {
    if (length(which(rownames(list@assays$RNA@counts) == tmp[j])) == 0){
      next() #Need to loop here because plot_hexbin_feature() does not have a built-in function to deal with absence of selected gene
    } else {
      plot <- plot_hexbin_feature(list, feature = tmp[j], type = "counts", action = "prop_0")+ 
        guides(fill=F, color = F) + 
        scale_fill_viridis(option = "B") 
      ggsave(path = "./UMAP/LineageMarkers", file = paste0(names(marker_list)[i], "_", tmp[j], "_prop.pdf"), plot, height=3, width=3.25)
    }
  }
}


#There are 20 gene sets I use as part of my manual annotation check (available in the
#./data/processedData directory as signature.list.rds). They are derived from
#[this](https://pubmed.ncbi.nlm.nih.gov/29961579/) publication. I use it as a baseline
#for my annotations and as a last slightly different approach to looking at cell type or
#possible phenotype using enrichment (this is different than using SingleR, ProjecTIL,
#or canonical markers). There are several methods for enrichment analysis, I am
#currently in the process of bringing them under one umbrella in the
#[escape](https://github.com/ncborcherding/escape) R package. However, I wanted to
#highlight work from Massimo Andreatta and Santiago Carmona - 
#[UCell](https://github.com/carmonalab/UCell) because with one core, it is faster than
#my original approach in escape and it is adapted to the issues with single-cell
#sequencing.




##############################################
#Loop Over the Enrichment Results
########################################

set.seed(123)
x <- sample(nrow(list[[]]), 20000)
set.names <- paste0(names(signature.list), "_UCell")
dir.create("./UMAP/GeneEnrichment")

for (i in seq_along(set.names)) {
  FeaturePlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 20000)], features =  set.names[i]) +
    scale_color_gradientn(colors = colorblind_vector(13)) + 
    theme(plot.title = element_blank())
  ggsave(paste0("./UMAP/GeneEnrichment/", set.names[i], ".pdf"), height = 3.5, width = 3.75)
}


# Contig Annotation

One issue with the data collection was the collation of TCR data from multiple sequencing runs of an [experiment](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469). I am including the code below for the sake of reproducibility, however, after isolating the individual TCR data, I manually deposited them into their respective directories.

data <- read.csv("./data/processedData/GSE144469_TCR_filtered_contig_annotations_all.csv")
data <- data[,-1]
data$sample <- stringr::str_split(data$barcode, "-", simplify = TRUE)[,2]
data$barcode<- paste0(stringr::str_split(data$barcode, "-", simplify = TRUE)[,1], "-1")

samples <- unique(data$sample)
for (i in seq_along(samples)) {
  tmp <- data[data$sample == samples[i],]
  tmp <- tmp[,-19] #remove the sample label that was added
  write.csv(tmp, paste0(samples[i], "_filtered_contig_annotations.csv"))
}


## Loading Contig Data and Sorting
library(scRepertoire)

######################################
#iterate to make a list of contig csvs
######################################

file_list <- list.files("./data/SequencingRuns", pattern = "filtered_contig_annotations", recursive = TRUE)
names <- stringr::str_split(file_list, "/", simplify = T)[,1]
contig.list <- NULL
for (i in seq_along(file_list)){
  contig.list[[i]] <-  read.csv(paste0("./data/SequencingRuns/", file_list[i]))
}
names(contig.list) <- names

##################################################
#Reducing the data to the individual barcode level
##################################################
combinedObject <- combineTCR(contig.list, samples = names, ID = rep("ID", length(contig.list)), filterMulti = TRUE, cells = "T-AB")

############################
#Refining the output a little
#############################

#Right now scRepertoire requires an ID variable to prevent issues with duplicate barcodes, this is not 
#ideal, but such is life, here I am just removing the ID to match the seurat object barcodes
#Yes I am the creator of scRepertoire, so this is really my own fault.
for (x in seq_along(combinedObject)) {
  combinedObject[[x]]$barcode <- stringr::str_remove_all(combinedObject[[x]]$barcode, "ID_")
}

names(combinedObject) <- stringr::str_remove_all(names(combinedObject), "_ID")

##############################################################
#Adding new variable "Patient" to calculate frequency by patient
################################################################
order.contigs <- names(combinedObject) #order of the names of the contigs
reorder.directory <- match(order.contigs, directory$SampleLabel) #matching names of order with directory
var <- directory$Sample[reorder.directory] #getting the patient order

combinedObject <- addVariable(combinedObject, name = "Patient", variables = var)


saveRDS(combinedObject, file = "./data/processedData/CombinedTCR_object.rds")

combinedObject <- readRDS("./data/processedData/CombinedTCR_object.rds")

list <- combineExpression(combinedObject, list, groupBy = "Patient") #Here is where the patient column comes in

list$cloneType <- factor(list$cloneType, levels = c("Rare (0 < X <= 1e-04)", "Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)", "Hyperexpanded (0.1 < X <= 1)"))

#Going to look at 80,000 cells in the UMAP just to better see distribution

DimPlot(list, cells = rownames(list[[]])[sample(nrow(list[[]]), 80000)], group.by  =  "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())
ggsave("./UMAP/ClonotypeExpansion.pdf", height = 3.5, width = 6)


saveRDS(list, file = "./data/ProcessedData/filtered_seuratObjects_harmony.rds")
