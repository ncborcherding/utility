# R/60_make_final_object.R

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(yaml)
  library(future)
  library(future.apply)
  library(scran)
  library(scRepertoire)
})

source("R/00_utils.R")

#' Generate the Full BPCells Directory
generate_full_bpcells <- function(config) {
  log_message("--- Generating FULL BPCells Dataset (Skipping Subset) ---")
  
  paths <- config$paths
  filter_config <- config$filtering
  
  # Output directory for FULL BPCells matrices
  bp_dir_full <- file.path(paths$results_dir, "01_bpcells_data_FULL")
  
  # If it exists, we assume it's done
  if (dir.exists(bp_dir_full)) {
    log_message("Full BPCells directory already exists. Using existing data.")
    return(bp_dir_full)
  }
  dir.create(bp_dir_full, recursive = TRUE)
  
  obj_dir <- paths$seurat_objects
  rds_files <- list.files(obj_dir, pattern = "\\.rds$", full.names = TRUE)
  
  for (file in rds_files) {
    sample_name <- tools::file_path_sans_ext(basename(file))
    
    obj <- readRDS(file)
    
    # Apply Filtering (Same as subset run)
    if (filter_config$run) {
      meta <- obj@meta.data
      cells_to_keep <- rep(TRUE, ncol(obj))
      for (f in filter_config$filters) {
        if (f$column %in% colnames(meta)) {
          col_data <- meta[[f$column]]
          if (f$type == "in") cells_to_keep <- cells_to_keep & (col_data %in% f$values)
          else if (f$type == "equals") cells_to_keep <- cells_to_keep & (col_data == f$values)
          else if (f$type == "not_na") cells_to_keep <- cells_to_keep & !is.na(col_data)
        }
      }
      if (sum(cells_to_keep) == 0) next
      obj <- subset(obj, cells = colnames(obj)[which(cells_to_keep)])
    }
    
    if (ncol(obj) == 0) next
    
    sample_bp_dir <- file.path(bp_dir_full, sample_name)
    
    counts_mat <- obj@assays$RNA@layers$counts
    if (is.null(counts_mat)) counts_mat <- obj@assays$RNA@counts
    
    BPCells::write_matrix_dir(mat = counts_mat, dir = sample_bp_dir)
    
    # Save metadata/features
    saveRDS(obj[[]], file = file.path(bp_dir_full, paste0(sample_name, "_meta.rds")))
    saveRDS(rownames(obj), file = file.path(bp_dir_full, paste0(sample_name, "_features.rds")))
    
    rm(obj, counts_mat); gc()
  }
  
  return(bp_dir_full)
}

make_final_object <- function(cfg, best_method, k, resolution) {
  
  log_message("=== Starting Creation of Final Full Object ===")
  log_message("Selected Best Method: ", best_method)
  log_message("Selected K: ", k, " | Resolution: ", resolution)
  
  # 1. Ensure Full Data Exists
  full_bp_dir <- generate_full_bpcells(cfg)
  
  # 2. Load Full Data
  log_message("Loading Full BPCells matrices...")
  sample_dirs <- list.dirs(full_bp_dir, full.names = FALSE, recursive = FALSE)
  sample_names <- sample_dirs[!grepl("_meta.rds|_features.rds", sample_dirs)]
  
  # Set up parallel plan for reading/processing
  n_workers <- cfg$memory$num_workers %||% 4
  if(n_workers == 0) n_workers <- future::availableCores() - 12
  plan("multisession", workers = n_workers)
  
  # Load objects into a list (lightweight pointers)
  object_list <- lapply(sample_names, function(s_name) {
    counts_bp <- BPCells::open_matrix_dir(file.path(full_bp_dir, s_name))
    meta <- readRDS(file.path(full_bp_dir, paste0(s_name, "_meta.rds")))
    rownames(counts_bp) <- readRDS(file.path(full_bp_dir, paste0(s_name, "_features.rds")))
    colnames(counts_bp) <- rownames(meta)
    Seurat::CreateSeuratObject(counts = counts_bp, meta.data = meta)
  })
  
  # 3. Preprocessing: Batch-aware HVG Calculation (scran style)
  log_message("Calculating HVGs per sample using scran::modelGeneVar (Parallel)...")
  
  # Use future_lapply for speed on 700+ samples
  per_sample_stats <- future_lapply(object_list, function(obj) {
    # Access BPCells counts
    counts_bp <- obj@assays$RNA@layers$counts
    if (is.null(counts_bp)) counts_bp <- obj@assays$RNA@counts
    
    # Calculate library sizes efficiently
    lib_sizes <- as.numeric(BPCells::colSums(counts_bp))
    lib_sizes[lib_sizes == 0] <- 1
    scale_factor <- median(lib_sizes)
    
    # Normalize (lazy op)
    norm_bp <- t(t(counts_bp) / lib_sizes) * scale_factor
    
    norm_mat <- as(norm_bp, "dgCMatrix")
    logcounts <- log1p(norm_mat)
    
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = logcounts))
    diff <- scran::modelGeneVar(sce, assay.type = "logcounts")
    diff@rownames <- rownames(obj)
    return(diff)
  }, future.seed = TRUE)
  
  # Combine variance stats
  log_message("Combining variance statistics...")
  common_genes <- Reduce(intersect, lapply(per_sample_stats, rownames))
  per_sample_stats <- lapply(per_sample_stats, function(stat) stat[common_genes, , drop = FALSE])
  combined_stats <- scran::combineVar(per_sample_stats)
  
  n_features <- cfg$preprocess$n_variable_features
  top_hvg <- scran::getTopHVGs(combined_stats, n = n_features)
  
  # Filter unwanted genes (MT, Ribosomal, VDJ)
  variable_features <- scRepertoire::quietVDJgenes(top_hvg)
  variable_features <- variable_features[!grepl('^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-', variable_features)]
  log_message("Identified ", length(variable_features), " robust variable features.")
  
  # 4. Merge and Scale
  log_message("Merging into single BPCells object...")
  if (length(object_list) > 1) {
    obj.full <- merge(object_list[[1]], y = object_list[2:length(object_list)])
  } else {
    obj.full <- object_list[[1]]
  }
  rm(object_list); gc()
  
  VariableFeatures(obj.full) <- variable_features
  
  log_message("Scaling and Running PCA on Full Dataset...")
  obj.full <- ScaleData(obj.full, features = variable_features)
  obj.full <- RunPCA(obj.full, npcs = 50, verbose = TRUE)
  
  # 5. Integration
  log_message("Running Integration: ", best_method)
  
  # Save preprocessed full object if needed by integration scripts
  full_preprocessed_path <- file.path(cfg$paths$results_dir, "03_preprocessed_data", "full_object_preprocessed.rds")
  safe_save_rds(obj.full, full_preprocessed_path)
  
  integration_func <- get(paste0("integrate_", best_method))
  full_integrated_path <- integration_func(full_preprocessed_path, "config.yaml")
  
  # Reload for clustering
  obj.integrated <- readRDS(full_integrated_path)
  
  # 6. Clustering (Optimized)
  log_message("Building Neighbor Graph (Annoy)...")
  
  reduction_use <- if(best_method == "harmony") "harmony" else 
    if(best_method %in% c("scvi", "scanvi")) "integrated.dr" else "pca"
  if(!reduction_use %in% names(obj.integrated@reductions)) reduction_use <- "pca"
  
  # Seurat FindNeighbors uses Annoy for efficiency
  obj.integrated <- FindNeighbors(
    obj.integrated, 
    reduction = reduction_use, 
    dims = 1:30, 
    k.param = k,
    return.neighbor = FALSE,
    method = "annoy" 
  )
  
  log_message("Running Leiden Clustering (Resolution ", resolution, ")...")
  obj.integrated <- FindClusters(
    obj.integrated, 
    resolution = resolution, 
    algorithm = 4
  )
  
  # 7. Save Final
  final_out_path <- file.path(cfg$paths$results_dir, "FINAL_integrated_object.rds")
  log_message("Saving Final Object to: ", final_out_path)
  safe_save_rds(obj.integrated, final_out_path)
  
  log_message("=== Final Object Created Successfully ===")
  return(final_out_path)
}