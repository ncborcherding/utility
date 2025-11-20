# R/10_process_seurat_runs.R
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(scDblFinder)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(scRepertoire)
})

source("R/00_utils.R")

qc_add_mito_ribo <- function(obj,
                             mito_pattern = "^MT-",
                             ribo_pattern = "^RPS|RPL-") {
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  obj
}

qc_filter <- function(obj, qc_cfg) {
  # Base guards
  keep <- rep(TRUE, ncol(obj))
  md   <- obj[[]]
  
  if (!is.null(qc_cfg$min_features)) {
    keep <- keep & md$nFeature_RNA >= qc_cfg$min_features
  }
  if (!is.null(qc_cfg$min_counts)) {
    keep <- keep & md$nCount_RNA   >= qc_cfg$min_counts
  }
  if (!is.null(md$percent.mito) && !is.null(qc_cfg$max_mito_pct)) {
    keep <- keep & md$percent.mito <= qc_cfg$max_mito_pct
  }
  if (!is.null(md$percent.ribo) && !is.null(qc_cfg$max_ribo_pct)) {
    keep <- keep & md$percent.ribo <= qc_cfg$max_ribo_pct
  }
  
  if (isTRUE(qc_cfg$high_feature_sd_enable %||% TRUE)) {
    sdmult <- qc_cfg$high_feature_sd %||% 2.5
    nf     <- md$nFeature_RNA
    if (length(nf) > 10 && all(nf > 0)) {
      cut <- round(exp(mean(log(nf)) + sdmult * sd(log(nf))))
      keep <- keep & (nf <= cut)
    }
  }
  
  subset(obj, cells = colnames(obj)[keep])
}

run_doublets <- function(obj, cores = 1) {
  bpp <- BiocParallel::MulticoreParam(workers = cores)
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, BPPARAM = bpp)
  df  <- data.frame(
    dbl.class = sce$scDblFinder.class,
    dbl.score = sce$scDblFinder.score,
    row.names = colnames(sce)
  )
  AddMetaData(obj, df)
}

process_seurat_runs <- function(config_path = "config.yaml") {
  
  cfg <- yaml::read_yaml(config_path)
  
  paths <- cfg$paths
  qc    <- cfg$qc
  dbl   <- cfg$doublets
  ann   <- cfg$annotation
  ss    <- cfg$subset
  
  results_dir <- paths$results_dir
  figures_dir <- paths$figures_dir
  tmp_dir     <- paths$tmp_dir
  
  tenx_dir <- paths$sequencing_runs
  
  tenx_runs <- if (!is.null(tenx_dir) && dir.exists(tenx_dir))
    list.dirs(tenx_dir, full.names = TRUE, recursive = FALSE) else character(0)
  
  for (src in tenx_runs) {
    m <- Read10X(src)
    if (inherits(m, "list")) m <- m[[1]]
    
    #Ensure using updated features
    updated.features <- UpdateGenes(rownames(m), gene.symbols = gene.symbols)
    # Sum counts for genes with the same (updated) symbol
    summed_matrix <-sparseRowsum(X = tmp, 
                                 group = updated.features[[1]], 
                                 drop.unused.levels = TRUE) 
    
    # Create a Seurat assay object from the summed count matrix
    tmp_assay <- CreateAssayObject(counts = summed_matrix)
    
    obj <- Createobject(counts = tmp_assay, assay = "RNA", project = basename(src))
    obj <- RenameCells(object = obj, 
                       new.names = paste0(basename(src), "_", rownames(obj[[]])))
    
    # Quality Control
    obj <- qc_add_mito_ribo(obj)
    obj <- qc_filter(obj, qc)
  
    # Estimating Doublets
    obj <- run_doublets(obj)
    
    # Annotating Cells
    Az <- Azimuth::RunAzimuth(obj, reference = "pbmcref")
    Az_annotation <- Az[[]][,grep("predicted", colnames(Az[[]]))]
    obj <- AddMetaData(obj, Az_annotation)
    
    sce   <- as.SingleCellExperiment(obj)
    HPCA  <- celldex::HumanPrimaryCellAtlasData()
    Monaco<- celldex::MonacoImmuneData()
    
    res1 <- SingleR(sce, ref = HPCA,   labels = HPCA$label.fine,   assay.type.test = 1)
    res2 <- SingleR(sce, ref = Monaco, labels = Monaco$label.fine, assay.type.test = 1)
    
    df1 <- data.frame(HPCA.labels = res1$labels,
                      HPCA.pruned = res1$pruned.labels,
                      row.names = rownames(res1))
    df2 <- data.frame(Monaco.labels = res2$labels,
                      Monaco.pruned = res2$pruned.labels,
                      row.names = rownames(res2))
    
    obj <- AddMetaData(obj, df1)
    obj <- AddMetaData(obj, df2)
    
    # Adding TCR
    TCR.file <- list.files(src, pattern = "annotations")[1]
    if(is.na(TCR.file)) {
      obj$CTaa <- NA
      obj$CTnt <- NA
      obj$CTgene <- NA
      obj$CTstrict <- NA
      obj$clonalProportion <- NA
      obj$clonalFrequency <- NA
      obj$cloneSize <- NA
    } else {
      TCR.file <- read.csv(file.path(src, TCR.file))
      combinedObject <- combineTCR(TCR.file, 
                                   samples = basename(src), 
                                   filterMulti = TRUE)
      obj <- combineExpression(combinedObject, 
                               obj, 
                               cloneCall = "aa")
    }
    
    saveRDS(obj, file = file.path(cfg$paths$seurat_objects, paste0(basename(src), ".rds")))
  }
  log_message("--- Individual Seurat Runs Processing Finished ---")
}