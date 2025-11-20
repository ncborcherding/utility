# R/35_integrate_scanvi.R
#
# Builds a minimal AnnData (cells x genes CSR counts) for scANVI,
# stamps Seurat HVGs into var['highly_variable'], and runs the Python step.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(yaml)
  library(reticulate)
})

source("R/00_utils.R")

integrate_scanvi <- function(preprocessed_obj_path, config_path = "config.yaml") {
  log_message("--- Starting scANVI Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)
  
  paths      <- config$paths
  labels_key <- config$methods$scanvi$labels_key
  batch_var  <- config$methods$harmony$group_by_vars[[1]] %||% NULL
  unlabeled  <- config$methods$scanvi$unlabeled_category %||% "Unknown"
  
  integration_timer <- start_timer() 
  
  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)
  DefaultAssay(obj) <- "RNA"
  set.seed(config$seed)
  
  if (!labels_key %in% colnames(obj@meta.data)) {
    stop("scANVI labels_key '", labels_key, "' not found in object metadata.")
  }
  
  # --- HVGs from Seurat (optional but recommended) ---
  hvg <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
  if (length(hvg) == 0L) message("No VariableFeatures() found; proceeding without precomputed HVGs.")
  
  # --- Build minimal obs (keep only label, batch, a few QC cols) ---
  keep_meta <- unique(c(
    labels_key,
    batch_var,
    intersect(c("nCount_RNA","nFeature_RNA","percent.mt"), colnames(obj@meta.data))
  ))
  obs_df <- as.data.frame(obj@meta.data[, keep_meta, drop = FALSE])
  # Ensure label column exists and coerce to character
  obs_df[[labels_key]] <- as.character(obs_df[[labels_key]])
  obs_df[[labels_key]][is.na(obs_df[[labels_key]]) | obs_df[[labels_key]] == ""] <- unlabeled
  if (!is.null(batch_var) && !(batch_var %in% names(obs_df))) {
    warning("Batch variable '", batch_var, "' not found; scANVI will run without batch_key.")
  }
  rownames(obs_df) <- colnames(obj)
  
  # --- Validate counts and names ---
  counts <- methods::as(SeuratObject::LayerData(obj), "dgCMatrix")
  if (is.null(counts) || nrow(counts) == 0)
    stop("RNA@counts is missing or empty. scANVI requires raw counts in X.")
  
  cells <- colnames(obj)
  genes <- rownames(obj)
  if (anyDuplicated(cells)) {
    stop("Duplicated cell names detected; make them unique first.")
  }
  genes_unique <- if (anyDuplicated(genes)) make.unique(genes) else genes
  
  # --- Reticulate env ---
  log_message("Using Conda env 'sc-integration-benchmark' via reticulate...")
  tryCatch({
    reticulate::use_condaenv("sc-integration-benchmark", required = TRUE)
    reticulate::py_config()
  }, error = function(e) {
    stop("Conda environment 'sc-integration-benchmark' not found. ", e$message)
  })
  
  # --- Write lean AnnData to tmp/scanvi_input.h5ad (cells x genes CSR counts) ---
  tmp_dir <- paths$tmp_dir
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  h5ad_path <- file.path(tmp_dir, "scanvi_input.h5ad")
  log_message("Writing AnnData (counts -> X) to: ", h5ad_path)
  
  anndata <- reticulate::import("anndata", convert = FALSE)
  sp      <- reticulate::import("scipy.sparse", convert = FALSE)
  np      <- reticulate::import("numpy", convert = FALSE)
  pd      <- reticulate::import("pandas", convert = TRUE)
  
  # dgCMatrix (CSC): p (indptr for columns), i (row indices), x (data)
  indptr_py  <- np$array(as.integer(counts@p), dtype = "int64")   # len = ncol+1 (cells+1)
  indices_py <- np$array(as.integer(counts@i), dtype = "int64")   # row indices (genes)
  data_py    <- np$array(as.numeric(counts@x), dtype = "float32")
  
  triplet     <- reticulate::tuple(data_py, indices_py, indptr_py)
  shape_tuple <- reticulate::tuple(as.integer(nrow(counts)), as.integer(ncol(counts)))  # genes x cells
  
  X_csc <- sp$csc_matrix(triplet, shape = shape_tuple)
  X     <- X_csc$transpose()$tocsr()  # cells x genes (CSR), as scvi/scanvi expect
  
  # obs/var
  obs_py <- pd$DataFrame(obs_df, index = cells)
  
  is_hvg <- genes_unique %in% hvg
  var_df <- data.frame(
    feature_name     = genes_unique,
    highly_variable  = is_hvg,
    stringsAsFactors = FALSE,
    check.names      = FALSE,
    row.names        = genes_unique
  )
  var_py <- pd$DataFrame(var_df)
  
  ad_obj <- anndata$AnnData(X = X, obs = obs_py, var = var_py)
  ad_obj$write_h5ad(h5ad_path, compression = "gzip")
  log_message("AnnData written successfully: ", h5ad_path)
  
  # --- Run Python scANVI step ---
  log_message("Sourcing Python scANVI script...")
  scanvi_py_script <- "py/scanvi_integration.py"
  if (!file.exists(scanvi_py_script)) {
    stop("Python script not found at: ", scanvi_py_script)
  }
  reticulate::source_python(scanvi_py_script)
  
  log_message("Training scANVI...")
  run_scanvi_integration(h5ad_path, config_path)  # modifies file in place
  
  # --- Load embeddings and attach to Seurat ---
  log_message("Loading scANVI embeddings from updated .h5ad ...")
  ad <- reticulate::import("anndata", convert = FALSE)
  ad_obj <- ad$read_h5ad(h5ad_path)
  
  scanvi_embed <- py_to_r(ad_obj$obsm["X_scanvi"])
  rownames(scanvi_embed) <- colnames(obj)
  
  obj[["scanvi"]] <- Seurat::CreateDimReducObject(
    embeddings = scanvi_embed,
    key = "scanvi_",
    assay = DefaultAssay(obj)
  )
  
  stop_timer(integration_timer, "scANVI Integration (save)")
  
  # Join Layers
  obj <- JoinLayers(obj)
  
  # --- Save ---
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "scanvi.rds")
  safe_save_rds(obj, output_path)
  
  log_message("Seurat object with scANVI integration saved: ", output_path)
  return(output_path)
}

# --- Script Execution guard (optional) ---
if (sys.nframe() == 0) {
  log_message("Running 35_integrate_scanvi.R as a standalone script.")
  log_message("Script 35_integrate_scanvi.R finished.")
}
