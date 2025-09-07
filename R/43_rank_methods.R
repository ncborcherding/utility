# R/40_rank_methods.R
# Aggregate all metrics across methods, rank them, and return the best method

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(yaml)
})

rank_methods <- function(config_path) {
  config <- yaml::read_yaml(config_path)
  
  # Load collected metrics from disk
  # (adjust these paths/names depending on how you save them)
  batch_metrics  <- readr::read_csv(file.path(config$paths$results_dir, "batch_metrics.csv"))
  label_metrics  <- readr::read_csv(file.path(config$paths$results_dir, "label_metrics.csv"))
  
  # Merge into a summary per method
  # Assumes each has a "method" column
  metrics <- dplyr::left_join(batch_metrics, label_metrics, by = "method")
  
  metrics <- metrics %>%
    mutate(
      z_batch = scale(batch_score)[,1],
      z_label = scale(label_score)[,1],
      score   = (z_batch + z_label) / 2
    )
  
  # Rank methods
  ranked <- metrics %>%
    arrange(desc(score))
  
  out_csv <- file.path(config$paths$results_dir, "method_ranking.csv")
  readr::write_csv(ranked, out_csv)
  
  # Log best method
  best_method <- ranked$method[1]
  message("Top-ranked integration method: ", best_method)
  message("Full ranking written to: ", out_csv)
  
  # Return best method for downstream clustering
  return(best_method)
}
