# run_pipeline.R
#
# Main orchestrator script for the integration benchmark pipeline.

# --- 1. Preamble ---

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# Source all R scripts in the R/ directory
# This makes their functions available to the pipeline runner.
r_scripts <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (script in r_scripts) {
  if (basename(script) != "install_packages.R" && basename(script) != "run_pipeline.R") {
    source(script)
  }
}
source("R/00_utils.R") # Source utils again to ensure it's available

# --- 2. Command-Line Argument Parsing ---

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yaml",
              help = "Path to the YAML configuration file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

config_path <- opts$config

# --- 3. Main Pipeline Function ---

run_pipeline <- function(config_path) {

  # --- Setup ---
  log_message("=== Starting Integration Benchmark Pipeline ===")

  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  config <- yaml::read_yaml(config_path)

  # Set seed for reproducibility
  set.seed(config$seed)

  # Create output directories
  dir.create(config$paths$results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(config$paths$figures_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(config$paths$tmp_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Pipeline Steps ---

  # Step 1: Process raw Seurat objects to on-disk BPCells format
  # This includes filtering and subsetting per-sample.
  bpcells_dir <- process_to_bpcells(config_path)

  # Step 2: Load BPCells data, find HVGs, merge, and preprocess (Scale/PCA)
  preprocessed_obj_path <- preprocess_bpcells_data(bpcells_dir, config_path)

  # Step 3: Run integration, analysis, and metrics for each method
  methods_to_run <- config$methods$run
  log_message("=== Starting Analysis for ", length(methods_to_run), " Methods ===")

  for (method in methods_to_run) {
    log_message(paste0(">>> Processing Method: ", toupper(method), " <<<"))

    # 4a. Integration
    integration_func <- get(paste0("integrate_", method))
    integrated_obj_path <- integration_func(preprocessed_obj_path, config_path)

    # 4b. Post-integration graph/UMAP
    analyzed_obj_path <- generate_graph_umap(integrated_obj_path, method, config_path)

    # 4c. Batch metrics
    calculate_batch_metrics(analyzed_obj_path, method, config_path)

    # 4d. Label metrics
    calculate_label_metrics(analyzed_obj_path, method, config_path)

    log_message(paste0(">>> Finished Processing: ", toupper(method), " <<<"))
  }

  # Step 5: Rank methods based on all collected metrics
  log_message("=== Aggregating Metrics and Ranking Methods ===")
  rank_methods(config_path)

  # Step 6: Generate final summary plots
  log_message("=== Generating Final Plots ===")
  generate_plots(config_path)

  log_message("=== Pipeline Finished Successfully ===")
}

# --- 4. Execute Pipeline ---

# Record start time
pipeline_timer <- start_timer()

# Run the main function
run_pipeline(config_path)

# Print total elapsed time
stop_timer(pipeline_timer, "Total pipeline execution")
