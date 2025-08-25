# R/40_graph_umap.R
#
# Computes graph, clusters, and UMAP for a given integration result.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
generate_graph_umap <- function(integrated_obj_path, method_name, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message(paste0("--- Generating Graph/UMAP for: ", method_name, " ---"))
  log_message("Reading configuration from: ", config_path)
  config <- yaml::read_yaml(config_path)

  params <- config$post_integration

  log_message("Loading integrated Seurat object from: ", integrated_obj_path)
  obj <- readRDS(integrated_obj_path)

  set.seed(config$seed)

  # --- 2. Check for Reduction ---
  # The reduction name should be the same as the method name (e.g., "harmony", "scvi")
  if (!method_name %in% Reductions(obj)) {
    stop("Reduction '", method_name, "' not found in the Seurat object at: ", integrated_obj_path)
  }

  log_message("Using reduction '", method_name, "' for downstream analysis.")

  # --- 3. Find Neighbors, Clusters, and UMAP ---
  analysis_timer <- start_timer()

  # a. Find Neighbors
  log_message("Finding neighbors...")
  obj <- FindNeighbors(
    obj,
    reduction = method_name,
    dims = 1:params$n_dims_use,
    k.param = params$n_neighbors,
    verbose = FALSE
  )

  # b. Find Clusters (Leiden)
  log_message("Finding clusters with Leiden algorithm...")
  # The name for the clusters metadata column will be specific to the method
  cluster_col_name <- paste0("leiden_", method_name)
  obj <- FindClusters(
    obj,
    resolution = params$leiden_resolution,
    algorithm = 4, # Leiden
    verbose = FALSE,
    graph.name = "RNA_snn" # Make sure to use the correct graph
  )
  # Rename the output column to be method-specific
  obj@meta.data[[cluster_col_name]] <- obj@meta.data$seurat_clusters
  obj@meta.data$seurat_clusters <- NULL # Remove the generic column

  # c. Run UMAP
  log_message("Running UMAP...")
  # The name for the UMAP reduction will be specific to the method
  umap_reduction_name <- paste0("umap_", method_name)
  obj <- RunUMAP(
    obj,
    reduction = method_name,
    dims = 1:params$n_dims_use,
    reduction.name = umap_reduction_name,
    verbose = FALSE
  )

  stop_timer(analysis_timer, paste0("Graph/UMAP generation for ", method_name))

  # --- 4. Save Updated Object and Clusters ---

  # Save the full object with the new data
  output_dir <- file.path(config$paths$results_dir, "05_analyzed_data")
  output_path_obj <- file.path(output_dir, paste0(method_name, "_analyzed.rds"))
  safe_save_rds(obj, output_path_obj)
  log_message("Analyzed Seurat object saved.")

  # Save the clusters to a separate CSV file
  output_path_csv <- file.path(output_dir, paste0("clusters_", method_name, ".csv"))
  clusters_df <- obj@meta.data %>%
    select(all_of(cluster_col_name)) %>%
    tibble::rownames_to_column("cell_id")

  write.csv(clusters_df, output_path_csv, row.names = FALSE)
  log_message("Clusters saved to CSV.")

  return(output_path_obj)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 40_graph_umap.R as a standalone script.")

  # This script is designed to be called by the main pipeline runner.
  # Example call:
  # generate_graph_umap("results/04_integrated_data/harmony.rds", "harmony", "config.yaml")

  log_message("Script 40_graph_umap.R finished.")
}
