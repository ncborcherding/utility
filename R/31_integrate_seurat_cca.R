# R/31_integrate_seurat_cca.R
#
# Runs the Seurat CCA integration workflow.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
integrate_seurat_cca <- function(preprocessed_obj_path, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message("--- Starting Seurat CCA Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)

  method_config <- config$methods$seurat_cca
  batch_var <- config$methods$harmony$group_by_vars[1] # Reuse harmony var as batch key

  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)

  set.seed(config$seed)

  # --- 2. Run Seurat CCA Integration ---
  integration_timer <- start_timer()

  # a. Split object by batch
  obj_list <- SplitObject(obj, split.by = batch_var)

  # b. Find Integration Anchors
  log_message("Finding integration anchors...")
  reduction_method <- ifelse(method_config$use_rpca, "rpca", "cca")
  log_message("Using reduction method: ", reduction_method)

  # Anchors are found on the raw 'SCT' or 'RNA' assay, using the pre-computed variable features
  # The default assay should be 'RNA' at this point.
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = VariableFeatures(obj),
    reduction = reduction_method,
    dims = 1:method_config$n_dims,
    verbose = FALSE
  )

  # c. Integrate Data
  log_message("Integrating data...")
  # This creates a new 'integrated' assay
  obj_integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:method_config$n_dims,
    verbose = FALSE
  )

  # --- 3. Post-Integration Processing ---
  log_message("Scaling integrated data and running PCA...")

  # Switch to the integrated assay for downstream analysis
  DefaultAssay(obj_integrated) <- "integrated"

  # Scale and run PCA on the integrated data
  obj_integrated <- ScaleData(obj_integrated, verbose = FALSE)
  obj_integrated <- RunPCA(
    obj_integrated,
    npcs = method_config$n_dims,
    verbose = FALSE
  )

  # To keep consistency with other methods that produce a single embedding,
  # we will store the final PCA results in a uniquely named DimReduc object
  # in the original object.

  # Rename the PCA reduction object
  obj_integrated@reductions$pca@key <- "cca_"
  colnames(obj_integrated@reductions$pca@cell.embeddings) <- paste0("cca_", 1:method_config$n_dims)

  # Copy the reduction back to the original object
  obj[["cca"]] <- obj_integrated@reductions$pca

  stop_timer(integration_timer, "Seurat CCA Integration")

  # --- 4. Save Integrated Object ---
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "seurat_cca.rds")

  safe_save_rds(obj, output_path)

  log_message("Seurat object with CCA integration saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 31_integrate_seurat_cca.R as a standalone script.")

  # Requires output from 30_preprocess.R
  # Example call:
  # integrate_seurat_cca("results/03_preprocessed_data/preprocessed_seurat_object.rds", "config.yaml")

  log_message("Script 31_integrate_seurat_cca.R finished.")
}
