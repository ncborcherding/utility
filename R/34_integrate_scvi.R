# R/34_integrate_scvi.R
#
# Runs the scVI integration workflow via reticulate.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(yaml)
  library(reticulate)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
integrate_scvi <- function(preprocessed_obj_path, config_path = "config.yaml") {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(yaml)
    library(reticulate)
    # SeuratDisk not needed anymore for conversion
  })
  
  source("R/00_utils.R")
  
  log_message("--- Starting scVI Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)
  
  paths <- config$paths
  # Keep your batch var lookup as-is; ensure it's present in meta later
  batch_var <- config$methods$harmony$group_by_vars[1]
  
  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)
  obj <- JoinLayers(obj)
  DefaultAssay(obj) <- "RNA"
  set.seed(config$seed)
  
  # ---- Validate counts and names ----
  counts <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
  if (is.null(counts) || nrow(counts) == 0)
    stop("RNA@counts is missing or empty. scVI requires raw counts in X.")
  
  if (!inherits(counts, "dgCMatrix"))
    counts <- as(counts, "dgCMatrix")
  
  cells <- colnames(obj)
  genes <- rownames(obj)
  
  if (anyDuplicated(cells)) {
    stop("Duplicated cell names detected in Seurat object. ",
         "Please make them unique before conversion (e.g., `Cells(obj) <- make.unique(Cells(obj))`).")
  }
  # AnnData dislikes duplicated var_names; make unique if needed (doesn't affect X order)
  genes_unique <- if (anyDuplicated(genes)) make.unique(genes) else genes
  
  # meta.data -> pandas.DataFrame (convert factors to strings for stability)
  obs_df <- as.data.frame(obj@meta.data)
  obs_df[] <- lapply(obs_df, function(x) if (is.factor(x)) as.character(x) else x)
  rownames(obs_df) <- cells
  if (!is.null(batch_var) && !(batch_var %in% colnames(obs_df))) {
    warning("Batch variable '", batch_var, "' not found in meta.data; ensure your Python step sets a valid batch key.")
  }
  
  log_message("Using Conda env 'sc-integration-benchmark' via reticulate...")
  tryCatch({
    reticulate::use_condaenv("sc-integration-benchmark", required = TRUE)
    reticulate::py_config()
  }, error = function(e) {
    stop("Conda environment 'sc-integration-benchmark' not found. ",
         "Create/activate it from your environment.yml. Error: ", e$message)
  })
  
  # ---- Build AnnData directly (proper tuples!) ----
  integration_timer <- start_timer()
  
  tmp_dir <- paths$tmp_dir
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  h5ad_path <- file.path(tmp_dir, "scvi_input.h5ad")
  log_message("Writing AnnData (counts -> X) to: ", h5ad_path)
  
  # ---- Build AnnData directly (proper orientation: cells x genes) ----
  anndata <- reticulate::import("anndata", convert = FALSE)
  sp      <- reticulate::import("scipy.sparse", convert = FALSE)
  np      <- reticulate::import("numpy", convert = FALSE)
  pd      <- reticulate::import("pandas", convert = TRUE)
  
  # dgCMatrix is CSC: p (indptr), i (indices), x (data)
  indptr_py  <- np$array(as.integer(counts@p), dtype = "int64")     # length ncol+1  (cells+1)
  indices_py <- np$array(as.integer(counts@i), dtype = "int64")     # row indices (genes)
  data_py    <- np$array(as.numeric(counts@x), dtype = "float32")
  
  # 3-tuple for (data, indices, indptr) and 2-tuple for shape
  triplet     <- reticulate::tuple(data_py, indices_py, indptr_py)
  shape_tuple <- reticulate::tuple(as.integer(nrow(counts)), as.integer(ncol(counts)))  # genes x cells
  
  # Build CSC (genes x cells), then transpose to CSR (cells x genes)
  X_csc <- sp$csc_matrix(triplet, shape = shape_tuple)
  X     <- X_csc$transpose()$tocsr()  # <-- key line; now rows = cells
  
  # Sanity checks (fail early if sizes don't line up)
  if (nrow(obj@meta.data) != ncol(counts)) {
    stop("obs rows (cells) = ", nrow(obj@meta.data), 
         " but counts has ", ncol(counts), " columns; they must match.")
  }
  if (length(genes_unique) != nrow(counts)) {
    stop("Number of genes in var (", length(genes_unique), 
         ") must equal nrow(counts) (", nrow(counts), ").")
  }
  
  # Build obs / var with indices matching X
  obs_df <- as.data.frame(obj@meta.data)
  obs_df[] <- lapply(obs_df, function(x) if (is.factor(x)) as.character(x) else x)
  rownames(obs_df) <- colnames(obj)  # cells
  obs_py <- pd$DataFrame(obs_df, index = colnames(obj))
  
  var_df <- data.frame(feature_name = genes_unique,
                       stringsAsFactors = FALSE, check.names = FALSE, row.names = genes_unique)
  var_py <- pd$DataFrame(var_df)
  
  ad_obj <- anndata$AnnData(X = X, obs = obs_py, var = var_py)
  ad_obj$write_h5ad(h5ad_path, compression = "gzip")
  log_message("AnnData written successfully.")

  # --- 3. Reticulate Setup and Python Script Execution ---
  log_message("Setting up reticulate to use Conda environment 'sc-integration-benchmark'...")
  tryCatch({
    use_condaenv("sc-integration-benchmark", required = TRUE)
    py_config() # Print config for debugging
  }, error = function(e) {
    stop("Conda environment 'sc-integration-benchmark' not found or reticulate failed. ",
         "Please ensure the environment is created and activated using 'environment.yml'. Error: ", e$message)
  })

  log_message("Executing Python script for scVI training...")

  # Pass arguments to the python script: h5ad path and config path
  # The python script will read its own parameters from the config file
  scvi_py_script <- "py/scvi_integration.py"
  if (!file.exists(scvi_py_script)) {
    stop("Python script not found at: ", scvi_py_script)
  }

  # Source the python script, which makes its functions available in R
  source_python(scvi_py_script)

  # Call the main function defined in the Python script
  # It will modify the .h5ad file in place
  run_scvi_integration(h5ad_path, config_path)

  log_message("Python script execution finished.")

  # --- 4. Load Results and Store Embeddings ---
  log_message("Loading scVI embeddings from updated .h5ad file...")

  # Load the AnnData object and extract the embeddings
  # Using anndata library directly via reticulate
  ad <- reticulate::import("anndata", convert = FALSE)
  ad_obj <- ad$read_h5ad(h5ad_path)

  # The python script will save embeddings in obsm
  scvi_embed <- py_to_r(ad_obj$obsm["X_scvi"])
  rownames(scvi_embed) <- rownames(obj@meta.data) # Ensure cell names are correct

  # Create a new DimReduc object
  scvi_reduction <- CreateDimReducObject(
    embeddings = scvi_embed,
    key = "scvi_",
    assay = DefaultAssay(obj)
  )

  # Add the new reduction to the Seurat object
  obj[["scvi"]] <- scvi_reduction

  stop_timer(integration_timer, "scVI Integration")

  # --- 5. Save Integrated Object ---
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "scvi.rds")

  safe_save_rds(obj, output_path)

  log_message("Seurat object with scVI integration saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 34_integrate_scvi.R as a standalone script.")

  # Requires output from 30_preprocess.R and the python script
  # Example call:
  # integrate_scvi("results/03_preprocessed_data/preprocessed_seurat_object.rds", "config.yaml")

  log_message("Script 34_integrate_scvi.R finished.")
}
