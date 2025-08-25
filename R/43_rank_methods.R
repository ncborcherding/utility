# R/43_rank_methods.R
#
# Ranks integration methods based on calculated metrics.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(yaml)
  library(tidyr)
  library(tibble)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
rank_methods <- function(config_path = "config.yaml") {

  # --- 1. Read Config and Load Metrics Data ---
  log_message("--- Ranking Integration Methods ---")
  config <- yaml::read_yaml(config_path)

  metrics_dir <- file.path(config$paths$results_dir, "06_metrics")
  if (!dir.exists(metrics_dir)) {
    stop("Metrics directory not found: ", metrics_dir)
  }

  rds_files <- list.files(metrics_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    stop("No metric files found in: ", metrics_dir)
  }

  log_message("Loading ", length(rds_files), " metric files.")
  all_metrics <- map_dfr(rds_files, readRDS)

  # --- 2. Normalize Metrics ---
  log_message("Normalizing metrics to a 0-1 scale (higher is better)...")

  # Pivot to a long format for easier processing
  metrics_long <- all_metrics %>%
    pivot_longer(
      cols = -method,
      names_to = "metric",
      values_to = "raw_value"
    )

  # Min-max normalization function
  min_max_norm <- function(x) {
    # If all values are the same, return 1 to avoid division by zero
    if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
      return(rep(1, length(x)))
    }
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }

  # Apply normalization per metric
  metrics_norm <- metrics_long %>%
    group_by(metric) %>%
    mutate(norm_value = min_max_norm(raw_value)) %>%
    ungroup()

  # --- 3. Compute Composite Scores ---
  log_message("Computing composite scores...")

  batch_metrics <- config$metrics$batch
  label_metrics <- config$metrics$label
  w_batch <- config$ranking$w_batch
  w_label <- config$ranking$w_label

  composite_scores <- metrics_norm %>%
    mutate(metric_type = case_when(
      metric %in% batch_metrics ~ "Batch",
      metric %in% label_metrics ~ "Label"
    )) %>%
    group_by(method, metric_type) %>%
    summarise(mean_norm_score = mean(norm_value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric_type,
      values_from = mean_norm_score,
      names_prefix = "Score_"
    ) %>%
    # Handle cases where a category of metrics wasn't run
    {
      if (!"Score_Batch" %in% names(.)) .$Score_Batch <- NA
      if (!"Score_Label" %in% names(.)) .$Score_Label <- NA
      .
    } %>%
    mutate(
      Score_Batch = coalesce(Score_Batch, 0),
      Score_Label = coalesce(Score_Label, 0),
      GlobalScore = (w_batch * Score_Batch) + (w_label * Score_Label)
    )

  # --- 4. Create Final League Table ---
  log_message("Creating final league table...")

  # Pivot normalized scores to wide format
  norm_wide <- metrics_norm %>%
    select(method, metric, norm_value) %>%
    pivot_wider(
      names_from = metric,
      values_from = norm_value,
      names_prefix = "norm_"
    )

  # Pivot raw scores to wide format
  raw_wide <- metrics_norm %>%
    select(method, metric, raw_value) %>%
    distinct() %>% # Avoid duplication from the long format
    pivot_wider(
      names_from = metric,
      values_from = raw_value
    )

  # Join everything together
  league_table <- raw_wide %>%
    left_join(norm_wide, by = "method") %>%
    left_join(composite_scores, by = "method") %>%
    arrange(desc(GlobalScore)) %>%
    mutate(Rank = row_number()) %>%
    # Reorder columns for clarity
    select(
      Rank,
      method,
      GlobalScore,
      Score_Batch,
      Score_Label,
      everything()
    )

  # --- 5. Save Table ---
  output_path <- file.path(config$paths$results_dir, "method_ranking.csv")
  log_message("Saving league table to: ", output_path)
  write.csv(league_table, output_path, row.names = FALSE)

  log_message("Method ranking finished.")
  return(output_path)
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 43_rank_methods.R as a standalone script.")
  # Example:
  # rank_methods("config.yaml")
  log_message("Script 43_rank_methods.R finished.")
}
