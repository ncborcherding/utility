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

# Ensure utils are loaded
source("R/00_utils.R") 
source("R/00_clust_utils.R") 

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
    #preprocessed_obj_path <- "results/03_preprocessed_data/preprocessed_bpcells_object.rds"
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
  best_method <- rank_methods(config_path)   
  best_method <- read.csv(best_method)[["method"]][1]
  log_message("Best integration method selected:", best_method)
  
  # Step 6: Generate final summary plots
  log_message("=== Generating Final Plots ===")
  generate_plots(config_path)
  
  log_message("=== Running Clustering Optimization on Best Method ===")
  
  # Locate integrated object from the best method
  best_integrated_obj_path <- file.path(
    config$paths$results_dir, "04_integrated_data",
    paste0(best_method, ".rds")
  )
  
  # Step 7: Grid search for leiden clustering
  run_cluster_param_grid(
    best_method = best_method,
    seurat_rds_path = best_integrated_obj_path,
    dims = seq(config$clustering$dims[[1]], config$clustering$dims[[2]]),
    graph_name = config$clustering$graph_name,
    resolutions = config$clustering$grid$resolutions,
    k_params = config$clustering$grid$k_params,
    n_workers = config$parallel$nworkers
  )
  
  # Step 8: Assessing leiden cluster stability
  run_cluster_stability(
    seurat_rds_path = best_integrated_obj_path,
    dims = seq(config$clustering$dims[[1]], config$clustering$dims[[2]]),
    graph_name = config$clustering$graph_name,
    subsample_frac = config$clustering$stability$subsample_frac,
    n_repeats = config$clustering$stability$n_repeats,
    n_workers = config$parallel$nworkers
  )
  
  # Step 9: Scoring leiden clusters
  combine_and_score(
    w_sil = config$clustering$scoring_weights$w_sil,
    w_mod = config$clustering$scoring_weights$w_mod,
    w_conn = config$clustering$scoring_weights$w_conn,
    w_stab = config$clustering$scoring_weights$w_stab,
    singleton_penalty = config$clustering$scoring_weights$singleton_penalty
  )
  
  apply_best_clusters
  
  
  
  #TODO Final Integration and clustering of full object issue #11 and 13 on github
  #apply_best_clusters(
  #  seurat_rds_in = ,
  #  dims = seq(config$clustering$dims[[1]], config$clustering$dims[[2]]),
  #  graph_name = config$clustering$graph_name
  #)
  
  #TODO Cell Annotation issue #12
  # - Canonical Marker Plots for T Cells
  # - Cell Type Annotation Plots
  # - Cluster Assignments
  
  #TODO Export to scanpy/scirpy issue #15 on github
  
  #TODO Run Cohort Summarization
  
 
  writeLines(capture.output(sessionInfo()), "/summary/sessionInfo.txt")
  
  log_message("=== Pipeline Finished Successfully ===")
}

# --- 4. Execute Pipeline ---

# Record start time
pipeline_timer <- start_timer()

# Run the main function
run_pipeline(config_path)

# Print total elapsed time
stop_timer(pipeline_timer, "Total pipeline execution")
