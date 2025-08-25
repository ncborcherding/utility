# R/30_preprocess.R
#
# Loads BPCells data, finds variable features, merges into a single
# memory-efficient Seurat object, and performs scaling and PCA.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(scran)
  library(scRepertoire) # For quietVDJgenes
  library(BiocParallel)
  library(BPCells)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
preprocess_bpcells_data <- function(bpcells_dir, config_path = "config.yaml") {

  # --- 1. Read Config and Find BPCells Dirs ---
  log_message("--- Preprocessing BPCells Data ---")
  config <- yaml::read_yaml(config_path)

  n_features <- config$preprocess$n_variable_features
  set.seed(config$seed)
  options(Seurat.object.assay.version = "v5")

  if (!dir.exists(bpcells_dir)) {
    stop("BPCells directory not found: ", bpcells_dir)
  }

  sample_dirs <- list.dirs(bpcells_dir, full.names = FALSE, recursive = FALSE)
  sample_names <- sample_dirs[!grepl("_meta.rds|_features.rds", sample_dirs)]

  if (length(sample_names) == 0) {
    stop("No sample directories found in BPCells path: ", bpcells_dir)
  }
  log_message("Found ", length(sample_names), " BPCells samples.")

  # --- 2. Create Object List and Find Variable Features ---
  preprocess_timer <- start_timer()

  # a. Create a list of Seurat objects from BPCells data
  log_message("Constructing list of BPCells-backed Seurat objects...")
  object_list <- lapply(sample_names, function(s_name) {
    counts_bp <- BPCells::open_matrix_dir(file.path(bpcells_dir, s_name))

    # Load metadata and features
    meta <- readRDS(file.path(bpcells_dir, paste0(s_name, "_meta.rds")))
    features <- readRDS(file.path(bpcells_dir, paste0(s_name, "_features.rds")))

    rownames(counts_bp) <- features
    colnames(counts_bp) <- rownames(meta)

    Seurat::CreateSeuratObject(counts = counts_bp, meta.data = meta)
  })
  names(object_list) <- sample_names

  # b. Run modelGeneVar on each object (scran workflow)
  log_message("Running scran::modelGeneVar on each sample...")
  if (requireNamespace("future", quietly = TRUE)) {
    plan("multisession", workers = max(1, future::availableCores() - 1))
  }

  object_list <- lapply(object_list, function(x) NormalizeData(x, verbose = FALSE))

  per_sample_stats <- lapply(object_list, function(x) {
    sce <- as.SingleCellExperiment(x, assay = "RNA")
    assayNames(sce)[assayNames(sce) == "data"] <- "logcounts"
    modelGeneVar(sce, assay.type = "logcounts")
  })

  # c. Combine variance and get top HVGs
  common_genes <- Reduce(intersect, lapply(per_sample_stats, rownames))
  per_sample_stats <- lapply(per_sample_stats, function(stat) stat[common_genes, , drop = FALSE])

  combined_stats <- combineVar(per_sample_stats)
  top_hvg <- getTopHVGs(combined_stats, n = n_features)

  # d. Filter unwanted genes
  variable_features <- scRepertoire::quietVDJgenes(top_hvg)
  variable_features <- variable_features[!grepl('^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-', variable_features)]
  log_message("Identified ", length(variable_features), " variable features.")

  # --- 3. Merge, Scale, and PCA ---
  log_message("Merging samples into a single BPCells-backed Seurat object...")

  # Merge the list of objects. This remains memory-efficient.
  if (length(object_list) > 1) {
    obj <- merge(object_list[[1]], y = object_list[2:length(object_list)])
  } else {
    obj <- object_list[[1]]
  }
  rm(object_list); gc()

  log_message("Merged object has ", ncol(obj), " cells.")

  VariableFeatures(obj) <- variable_features

  log_message("Scaling data and running PCA...")
  obj <- ScaleData(obj, features = variable_features, verbose = FALSE)
  obj <- RunPCA(obj, features = variable_features, verbose = FALSE, npcs = 50)

  stop_timer(preprocess_timer, "BPCells Preprocessing")

  # --- 4. Save Preprocessed Object ---
  output_dir <- file.path(config$paths$results_dir, "03_preprocessed_data")
  output_path <- file.path(output_dir, "preprocessed_bpcells_object.rds")

  safe_save_rds(obj, output_path)
  log_message("Preprocessed BPCells-backed Seurat object saved.")

  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 30_preprocess.R as a standalone script.")
  # Example call:
  # preprocess_bpcells_data("results/01_bpcells_data", "config.yaml")
  log_message("Script 30_preprocess.R finished.")
}
