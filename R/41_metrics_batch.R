# R/41_metrics_batch.R
#
# Computes batch-mixing metrics for a given integration result.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(lisi)
  library(cluster)
  # For kBET, we can use a helper function to avoid installing the full kBET package
  # which can have heavy dependencies. Let's use the FNN package for kNN search.
  library(FNN)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- kBET Helper Function ---
# A simplified implementation of the kBET concept.
# Calculates the observed vs. expected batch frequencies in local neighborhoods.
calculate_kbet <- function(emb, batch, k = 25) {
  if (is.factor(batch)) batch <- droplevels(batch)

  # 1. Find k-nearest neighbors for each cell
  knn_obj <- FNN::get.knn(emb, k = k)

  # 2. For each cell's neighborhood, get the batch labels
  batch_freq_observed <- apply(knn_obj$nn.index, 1, function(indices) {
    table(batch[indices])
  })

  # 3. Calculate global batch frequencies
  global_freq <- table(batch) / length(batch)

  # 4. Perform chi-squared test for each neighborhood
  # This is a simplification. The official kBET does more complex things.
  # Here we just calculate the rejection rate of a simple test.
  # A higher p-value means we cannot reject the null (it's well-mixed).

  results <- sapply(1:nrow(emb), function(i) {
    # Get observed counts for the neighborhood of cell i
    obs_counts <- rep(0, length(levels(batch)))
    names(obs_counts) <- levels(batch)

    neighborhood_indices <- knn_obj$nn.index[i, ]
    neighborhood_batches <- batch[neighborhood_indices]

    tbl <- table(neighborhood_batches)
    obs_counts[names(tbl)] <- tbl

    # Expected counts based on global frequencies
    expected_counts <- global_freq * k

    # Chi-squared test
    # Add a small epsilon to avoid division by zero or errors with 0 counts
    test_result <- suppressWarnings(
      chisq.test(x = obs_counts + 1e-6, p = global_freq + 1e-6)
    )

    return(test_result$p.value)
  })

  # Rejection rate at alpha = 0.05
  rejection_rate <- mean(results < 0.05, na.rm = TRUE)

  # Return the "acceptance rate" so higher is better
  return(1 - rejection_rate)
}


# --- Main function ---
calculate_batch_metrics <- function(analyzed_obj_path, method_name, config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message(paste0("--- Calculating Batch Metrics for: ", method_name, " ---"))
  config <- yaml::read_yaml(config_path)

  metrics_to_run <- config$metrics$batch
  batch_var <- config$methods$harmony$group_by_vars[1] # Main batch key
  n_dims <- config$post_integration$n_dims_use

  obj <- readRDS(analyzed_obj_path)

  # --- 2. Extract Data ---
  log_message("Extracting embedding and metadata...")
  if (!method_name %in% Reductions(obj)) {
    stop("Reduction '", method_name, "' not found in the object.")
  }
  emb <- Embeddings(obj, reduction = method_name)[, 1:n_dims]

  meta <- obj@meta.data
  if (!batch_var %in% colnames(meta)) {
    stop("Batch variable '", batch_var, "' not found in metadata.")
  }
  batch_labels <- meta[[batch_var]]

  # --- 3. Calculate Metrics ---
  metrics_timer <- start_timer()
  results <- list()

  # ASW (Silhouette)
  if ("asw_batch" %in% metrics_to_run) {
    log_message("Calculating ASW (batch)...")
    # Subsample for speed if dataset is large
    n_cells <- nrow(emb)
    if (n_cells > 5000) {
      sample_indices <- sample(1:n_cells, 5000)
      emb_sample <- emb[sample_indices, ]
      batch_labels_sample <- batch_labels[sample_indices]
    } else {
      emb_sample <- emb
      batch_labels_sample <- batch_labels
    }

    dist_matrix <- dist(emb_sample)
    sil <- silhouette(as.integer(as.factor(batch_labels_sample)), dist_matrix)
    asw_score <- mean(sil[, "sil_width"])
    results$asw_batch <- 1 - asw_score # Higher is better mixing
  }

  # iLISI
  if ("ilisi" %in% metrics_to_run) {
    log_message("Calculating iLISI...")
    ilisi_scores <- lisi::compute_lisi(emb, meta, batch_var)
    results$ilisi <- mean(ilisi_scores[, 1], na.rm = TRUE)
  }

  # kBET
  if ("kbet_batch" %in% metrics_to_run) {
    log_message("Calculating kBET (batch)...")
    # kBET can be slow, subsample to 5k cells if larger
    n_cells <- nrow(emb)
    if (n_cells > 5000) {
      sample_indices <- sample(1:n_cells, 5000)
      emb_sample <- emb[sample_indices, ]
      batch_labels_sample <- as.factor(batch_labels[sample_indices])
    } else {
      emb_sample <- emb
      batch_labels_sample <- as.factor(batch_labels)
    }
    results$kbet_batch <- calculate_kbet(emb_sample, batch_labels_sample)
  }

  stop_timer(metrics_timer, paste0("Batch metrics calculation for ", method_name))

  # --- 4. Save Results ---
  metrics_df <- as.data.frame(results)
  metrics_df$method <- method_name

  output_dir <- file.path(config$paths$results_dir, "06_metrics")
  output_path <- file.path(output_dir, paste0("metrics_batch_", method_name, ".rds"))

  safe_save_rds(metrics_df, output_path)
  log_message("Batch metrics saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 41_metrics_batch.R as a standalone script.")
  # Example:
  # calculate_batch_metrics("results/05_analyzed_data/harmony_analyzed.rds", "harmony", "config.yaml")
  log_message("Script 41_metrics_batch.R finished.")
}
