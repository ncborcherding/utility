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

  # --- 1. Read Config and Load Data ---
  log_message("--- Starting scVI Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)

  paths <- config$paths
  batch_var <- config$methods$harmony$group_by_vars[1]

  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)

  set.seed(config$seed)

  # --- 2. Data Conversion (Seurat -> AnnData) ---
  integration_timer <- start_timer()

  # Create a temporary directory for h5ad file if it doesn't exist
  tmp_dir <- paths$tmp_dir
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  h5ad_path <- file.path(tmp_dir, "scvi_input.h5ad")
  log_message("Converting Seurat object to AnnData file at: ", h5ad_path)

  # scVI works on raw counts. We need to ensure the 'counts' slot is present.
  # The object from preprocessing should still have its original counts.
  DefaultAssay(obj) <- "RNA"
  SaveH5AD(obj, h5ad_path, assay = "RNA", overwrite = TRUE)
  log_message("Conversion complete.")

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
