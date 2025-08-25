# R/32_integrate_harmony.R
#
# Runs the Harmony integration workflow.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(harmony)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
integrate_harmony <- function(preprocessed_obj_path, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message("--- Starting Harmony Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)

  method_config <- config$methods$harmony

  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)

  set.seed(config$seed)

  # --- 2. Run Harmony Integration ---
  integration_timer <- start_timer()

  # Harmony runs on a PCA reduction. The preprocessed object should have this.
  if (!"pca" %in% Reductions(obj)) {
    stop("PCA reduction not found in the preprocessed object. Please run PCA first.")
  }

  # Get PCA embeddings and metadata
  pca_embed <- Embeddings(obj, reduction = "pca")
  meta_data <- obj@meta.data

  # Check if batch variables exist
  if (!all(method_config$group_by_vars %in% colnames(meta_data))) {
    stop("Not all Harmony 'group_by_vars' found in object metadata.")
  }

  log_message("Running Harmony on ", ncol(pca_embed), " PCs...")

  # Use the native HarmonyMatrix function
  harmony_embed <- HarmonyMatrix(
    data_mat  = pca_embed,
    meta_data = meta_data,
    vars_use  = method_config$group_by_vars,
    theta     = method_config$theta,
    lambda    = method_config$lambda,
    max.iter.harmony = method_config$max_iter_harmony,
    do_pca    = FALSE, # We are passing in a PCA, so no need to recalculate
    verbose   = FALSE
  )

  # --- 3. Store Harmony Embeddings ---
  log_message("Storing Harmony embeddings in the Seurat object...")

  # Create a new DimReduc object for the Harmony results
  harmony_reduction <- CreateDimReducObject(
    embeddings = harmony_embed,
    key = "harmony_",
    assay = DefaultAssay(obj)
  )

  # Add the new reduction to the Seurat object
  obj[["harmony"]] <- harmony_reduction

  stop_timer(integration_timer, "Harmony Integration")

  # --- 4. Save Integrated Object ---
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "harmony.rds")

  safe_save_rds(obj, output_path)

  log_message("Seurat object with Harmony integration saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 32_integrate_harmony.R as a standalone script.")

  # Requires output from 30_preprocess.R
  # Example call:
  # integrate_harmony("results/03_preprocessed_data/preprocessed_seurat_object.rds", "config.yaml")

  log_message("Script 32_integrate_harmony.R finished.")
}
