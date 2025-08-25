# R/42_metrics_label.R
#
# Computes cell label conservation metrics for a given integration result.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(lisi)
  library(aricode) # For NMI
  library(mclust)  # For ARI
  library(FNN)
})

# --- Source utilities ---
source("R/00_utils.R")
# Source the kBET helper from the other metrics script
source("R/41_metrics_batch.R") # This contains calculate_kbet

# --- Main function ---
calculate_label_metrics <- function(analyzed_obj_path, method_name, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message(paste0("--- Calculating Label Metrics for: ", method_name, " ---"))
  config <- yaml::read_yaml(config_path)

  metrics_to_run <- config$metrics$label
  labels_key <- config$methods$scanvi$labels_key
  n_dims <- config$post_integration$n_dims_use

  obj <- readRDS(analyzed_obj_path)

  # --- 2. Extract Data ---
  log_message("Extracting embedding and metadata...")
  if (!method_name %in% Reductions(obj)) {
    stop("Reduction '", method_name, "' not found in the object.")
  }
  emb <- Embeddings(obj, reduction = method_name)[, 1:n_dims]

  meta <- obj@meta.data
  cluster_col <- paste0("leiden_", method_name)

  # Check for required columns
  if (!labels_key %in% colnames(meta)) {
    stop("Ground truth label column '", labels_key, "' not found in metadata.")
  }
  if (!cluster_col %in% colnames(meta)) {
    stop("Cluster column '", cluster_col, "' not found for method '", method_name, "'.")
  }

  true_labels <- meta[[labels_key]]
  cluster_labels <- meta[[cluster_col]]

  # --- 3. Calculate Metrics ---
  metrics_timer <- start_timer()
  results <- list()

  # NMI (Normalized Mutual Information)
  if ("nmi" %in% metrics_to_run) {
    log_message("Calculating NMI...")
    results$nmi <- aricode::NMI(cluster_labels, true_labels)
  }

  # ARI (Adjusted Rand Index)
  if ("ari" %in% metrics_to_run) {
    log_message("Calculating ARI...")
    results$ari <- mclust::adjustedRandIndex(cluster_labels, true_labels)
  }

  # cLISI (cell-type LISI)
  if ("clisi" %in% metrics_to_run) {
    log_message("Calculating cLISI...")
    clisi_scores <- lisi::compute_lisi(emb, meta, labels_key)
    # Inverse score so higher = better (less mixing of labels)
    # Add a small epsilon to avoid division by zero if mean is 1
    mean_clisi <- mean(clisi_scores[, 1], na.rm = TRUE)
    results$clisi <- 1 / (mean_clisi + 1e-9)
  }

  # kBET (per label)
  if ("kbet_label" %in% metrics_to_run) {
    log_message("Calculating kBET (label)...")
    # Subsample for speed
    n_cells <- nrow(emb)
    if (n_cells > 5000) {
      sample_indices <- sample(1:n_cells, 5000)
      emb_sample <- emb[sample_indices, ]
      labels_sample <- as.factor(true_labels[sample_indices])
    } else {
      emb_sample <- emb
      labels_sample <- as.factor(true_labels)
    }

    # High acceptance rate means neighborhoods are pure by label, which is good.
    results$kbet_label <- calculate_kbet(emb_sample, labels_sample)
  }

  stop_timer(metrics_timer, paste0("Label metrics calculation for ", method_name))

  # --- 4. Save Results ---
  metrics_df <- as.data.frame(results)
  metrics_df$method <- method_name

  output_dir <- file.path(config$paths$results_dir, "06_metrics")
  output_path <- file.path(output_dir, paste0("metrics_label_", method_name, ".rds"))

  safe_save_rds(metrics_df, output_path)
  log_message("Label metrics saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 42_metrics_label.R as a standalone script.")
  # Example:
  # calculate_label_metrics("results/05_analyzed_data/harmony_analyzed.rds", "harmony", "config.yaml")
  log_message("Script 42_metrics_label.R finished.")
}
