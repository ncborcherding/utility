# R/33_integrate_fastmnn.R
#
# Runs the fastMNN integration workflow from the 'batchelor' package.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(batchelor)
  library(SingleCellExperiment)
  library(BiocParallel)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
integrate_fastmnn <- function(preprocessed_obj_path, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message("--- Starting fastMNN Integration ---")
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)

  method_config <- config$methods$fastmnn

  log_message("Loading preprocessed Seurat object from: ", preprocessed_obj_path)
  obj <- readRDS(preprocessed_obj_path)
  obj <- JoinLayers(obj) # Join layers to extract for as.SingleCellExperiment()

  set.seed(config$seed)

  # --- 2. Convert to SingleCellExperiment and Run fastMNN ---
  integration_timer <- start_timer()

  # fastMNN works on log-normalized counts of variable features.
  # Convert Seurat object to SCE
  log_message("Converting Seurat object to SingleCellExperiment...")

  # Use the pre-selected variable features
  hvg_list <- VariableFeatures(obj)
  if (length(hvg_list) == 0) {
    stop("Variable features not found in the preprocessed object.")
  }

  # Subset the object to only include variable features for the conversion
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  assayNames(sce)[assayNames(sce) == "data"] <- "logcounts" # Ensure logcounts are named correctly
  sce <- sce[hvg_list, ]

  # Check if batch variable exists
  batch_var <- method_config$batch_var
  if (!batch_var %in% colnames(colData(sce))) {
    stop("fastMNN 'batch_var' (", batch_var, ") not found in object metadata.")
  }

  # Set up parallel processing
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    n_workers <- max(1, future::availableCores() - 1)
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers)
    log_message("Using BiocParallel with ", n_workers, " workers.")
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }

  log_message("Running fastMNN...")

  # Run fastMNN on the SCE object, using the batch variable
  # This performs cosine normalization and multi-batch PCA internally
  sce_corrected <- fastMNN(
    sce,
    batch = colData(sce)[[batch_var]],
    d = method_config$n_dims,
    BPPARAM = BPPARAM
  )

  # --- 3. Store fastMNN Embeddings ---
  log_message("Storing fastMNN embeddings in the Seurat object...")

  # The corrected embeddings are in the 'corrected' reducedDim slot
  fastmnn_embed <- reducedDim(sce_corrected, "corrected")
  rownames(fastmnn_embed) <- colnames(sce_corrected)

  # Create a new DimReduc object for the fastMNN results
  fastmnn_reduction <- CreateDimReducObject(
    embeddings = fastmnn_embed,
    key = "fastmnn_",
    assay = DefaultAssay(obj)
  )

  # Add the new reduction to the original Seurat object
  obj[["fastmnn"]] <- fastmnn_reduction

  stop_timer(integration_timer, "fastMNN Integration")

  # --- 4. Save Integrated Object ---
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "fastmnn.rds")

  safe_save_rds(obj, output_path)

  log_message("Seurat object with fastMNN integration saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 33_integrate_fastmnn.R as a standalone script.")

  # Requires output from 30_preprocess.R
  # Example call:
  # integrate_fastmnn("results/03_preprocessed_data/preprocessed_seurat_object.rds", "config.yaml")

  log_message("Script 33_integrate_fastmnn.R finished.")
}
