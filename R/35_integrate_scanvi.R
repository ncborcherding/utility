# R/35_integrate_scanvi.R
#
# Runs the scANVI integration workflow via reticulate.

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
integrate_scanvi <- function(preprocessed_obj_path, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message("--- Starting scANVI Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)

  paths <- config$paths
  labels_key <- config$methods$scanvi$labels_key

  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)

  # Ensure the label column exists
  if (!labels_key %in% colnames(obj@meta.data)) {
    stop("scANVI labels_key '", labels_key, "' not found in object metadata.")
  }

  set.seed(config$seed)

  # --- 2. Data Conversion (Seurat -> AnnData) ---
  integration_timer <- start_timer()

  tmp_dir <- paths$tmp_dir
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  h5ad_path <- file.path(tmp_dir, "scanvi_input.h5ad")
  log_message("Converting Seurat object to AnnData file at: ", h5ad_path)

  DefaultAssay(obj) <- "RNA"
  SaveH5AD(obj, h5ad_path, assay = "RNA", overwrite = TRUE)
  log_message("Conversion complete.")

  # --- 3. Reticulate Setup and Python Script Execution ---
  log_message("Setting up reticulate to use Conda environment 'sc-integration-benchmark'...")
  tryCatch({
    use_condaenv("sc-integration-benchmark", required = TRUE)
    py_config()
  }, error = function(e) {
    stop("Conda environment 'sc-integration-benchmark' not found or reticulate failed. ",
         "Please ensure the environment is created and activated using 'environment.yml'. Error: ", e$message)
  })

  log_message("Executing Python script for scANVI training...")

  scanvi_py_script <- "py/scanvi_integration.py"
  if (!file.exists(scanvi_py_script)) {
    stop("Python script not found at: ", scanvi_py_script)
  }

  source_python(scanvi_py_script)

  # Call the main function from the Python script
  run_scanvi_integration(h5ad_path, config_path)

  log_message("Python script execution finished.")

  # --- 4. Load Results and Store Embeddings ---
  log_message("Loading scANVI embeddings from updated .h5ad file...")

  ad <- reticulate::import("anndata", convert = FALSE)
  ad_obj <- ad$read_h5ad(h5ad_path)

  # The python script will save embeddings in obsm with key "X_scanvi"
  scanvi_embed <- py_to_r(ad_obj$obsm["X_scanvi"])
  rownames(scanvi_embed) <- colnames(obj)

  # Create a new DimReduc object
  scanvi_reduction <- CreateDimReducObject(
    embeddings = scanvi_embed,
    key = "scanvi_",
    assay = DefaultAssay(obj)
  )

  # Add the new reduction to the Seurat object
  obj[["scanvi"]] <- scanvi_reduction

  stop_timer(integration_timer, "scANVI Integration")

  # --- 5. Save Integrated Object ---
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "scanvi.rds")

  safe_save_rds(obj, output_path)

  log_message("Seurat object with scANVI integration saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 35_integrate_scanvi.R as a standalone script.")

  # Requires output from 30_preprocess.R and the python script
  # Example call:
  # integrate_scanvi("results/03_preprocessed_data/preprocessed_seurat_object.rds", "config.yaml")

  log_message("Script 35_integrate_scanvi.R finished.")
}
