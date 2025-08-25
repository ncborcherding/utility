# R/15_process_to_bpcells.R
#
# Processes individual Seurat objects: filters, subsets, and converts to on-disk BPCells format.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(BPCells)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Helper function for filtering ---
apply_filters <- function(obj, filter_list) {
  meta <- obj@meta.data
  cells_to_keep <- rep(TRUE, ncol(obj))

  for (f in filter_list) {
    if (!f$column %in% colnames(meta)) {
      log_message("Warning: Filter column '", f$column, "' not found. Skipping this filter.")
      next
    }

    col_data <- meta[[f$column]]

    if (f$type == "in") {
      cells_to_keep <- cells_to_keep & (col_data %in% f$values)
    } else if (f$type == "equals") {
      cells_to_keep <- cells_to_keep & (col_data == f$values)
    } else if (f$type == "not_na") {
      cells_to_keep <- cells_to_keep & !is.na(col_data)
    }
  }

  return(subset(obj, cells = colnames(obj)[which(cells_to_keep)]))
}


# --- Main function ---
process_to_bpcells <- function(config_path = "config.yaml") {

  # --- 1. Read Config and Find Files ---
  log_message("--- Processing Seurat Objects to BPCells Format ---")
  config <- yaml::read_yaml(config_path)

  paths <- config$paths
  filter_config <- config$filtering
  subset_config <- config$subset
  set.seed(config$seed)

  # Input directory for raw Seurat objects
  obj_dir <- paths$seurat_objects
  if (!dir.exists(obj_dir)) stop("Seurat object directory not found: ", obj_dir)
  rds_files <- list.files(obj_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) stop("No RDS files found in: ", obj_dir)

  # Output directory for BPCells matrices
  bp_dir <- file.path(paths$results_dir, "01_bpcells_data")
  if (dir.exists(bp_dir)) {
    log_message("BPCells directory already exists. Deleting to ensure a fresh run.")
    unlink(bp_dir, recursive = TRUE)
  }
  dir.create(bp_dir, recursive = TRUE)

  # --- 2. Loop Through Objects for Filtering and Subsetting ---
  log_message("Starting loop over ", length(rds_files), " files...")

  # If subsetting, first get total cell counts after filtering to determine proportions
  total_cells_after_filter <- 0
  file_cell_counts <- list()
  if (subset_config$run) {
    log_message("Pre-scanning files to determine subsetting proportions...")
    for (file in rds_files) {
      obj <- readRDS(file)
      if (filter_config$run) {
        obj <- apply_filters(obj, filter_config$filters)
      }
      file_cell_counts[[basename(file)]] <- ncol(obj)
      total_cells_after_filter <- total_cells_after_filter + ncol(obj)
    }
  }

  for (file in rds_files) {
    sample_name <- tools::file_path_sans_ext(basename(file))
    log_message("Processing sample: ", sample_name)

    obj <- readRDS(file)

    # a. Apply filtering
    if (filter_config$run) {
      n_before <- ncol(obj)
      obj <- apply_filters(obj, filter_config$filters)
      log_message("Filtered from ", n_before, " to ", ncol(obj), " cells.")
    }

    if (ncol(obj) == 0) {
      log_message("Skipping sample ", sample_name, " as it has no cells after filtering.")
      next
    }

    # b. Apply subsetting
    if (subset_config$run) {
      n_before_subset <- ncol(obj)
      prop_cells <- file_cell_counts[[basename(file)]] / total_cells_after_filter
      n_to_sample <- round(prop_cells * subset_config$n_cells_total)

      if (n_to_sample < n_before_subset) {
        cells_to_keep <- sample(colnames(obj), n_to_sample)
        obj <- subset(obj, cells = cells_to_keep)
        log_message("Subsetted from ", n_before_subset, " to ", ncol(obj), " cells.")
      }
    }

    if (ncol(obj) == 0) {
      log_message("Skipping sample ", sample_name, " as it has no cells after subsetting.")
      next
    }

    # c. Write to BPCells format
    sample_bp_dir <- file.path(bp_dir, sample_name)

    # Using the 'counts' layer from the 'RNA' assay
    counts_mat <- obj@assays$RNA@layers$counts
    if (is.null(counts_mat)) {
        counts_mat <- obj@assays$RNA@counts
    }

    BPCells::write_matrix_dir(mat = counts_mat, dir = sample_bp_dir)

    # Save metadata and feature names
    saveRDS(obj[[]], file = file.path(bp_dir, paste0(sample_name, "_meta.rds")))
    saveRDS(rownames(obj), file = file.path(bp_dir, paste0(sample_name, "_features.rds")))

    rm(obj, counts_mat); gc()
  }

  log_message("--- BPCells Processing Finished ---")
  return(bp_dir)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 15_process_to_bpcells.R as a standalone script.")
  process_to_bpcells("config.yaml")
  log_message("Script 15_process_to_bpcells.R finished.")
}
