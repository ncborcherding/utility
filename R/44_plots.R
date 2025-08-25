# R/44_plots.R
#
# Generates all summary plots for the integration benchmark.

# --- Load libraries ---
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(purrr)
  library(tidyr)
  library(patchwork)
  library(gt)
  library(fmsb) # For radar charts
  library(RColorBrewer)
})

# --- Source utilities ---
source("R/00_utils.R")

# --- Main function ---
generate_plots <- function(config_path = "config.yaml") {

  # --- 1. Read Config and Load Data ---
  log_message("--- Generating Summary Plots ---")
  config <- yaml::read_yaml(config_path)

  paths <- config$paths
  ranking_file <- file.path(paths$results_dir, "method_ranking.csv")
  analyzed_dir <- file.path(paths$results_dir, "05_analyzed_data")

  if (!file.exists(ranking_file)) {
    stop("Ranking file not found: ", ranking_file)
  }

  league_table <- read.csv(ranking_file)

  # --- 2. Bar Plot of Raw Metrics ---
  log_message("Generating bar plot of raw metrics...")

  raw_metrics_long <- league_table %>%
    select(method, all_of(config$metrics$batch), all_of(config$metrics$label)) %>%
    pivot_longer(-method, names_to = "metric", values_to = "score")

  p_bar <- ggplot(raw_metrics_long, aes(x = reorder(method, -score), y = score, fill = method)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Raw Metric Scores by Integration Method",
      x = "Method",
      y = "Raw Score (directionally adjusted)"
    ) +
    UtilityTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(paths$figures_dir, "01_barplot_raw_metrics.png"), p_bar, width = 10, height = 8, bg = "white")

  # --- 3. Heatmap of Normalized Metrics ---
  log_message("Generating heatmap of normalized metrics...")

  norm_metrics_long <- league_table %>%
    select(method, starts_with("norm_")) %>%
    pivot_longer(-method, names_to = "metric", values_to = "score") %>%
    mutate(metric = gsub("norm_", "", metric))

  p_heatmap <- ggplot(norm_metrics_long, aes(x = metric, y = method, fill = score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(score, 2)), color = "black", size = 3) +
    scale_fill_viridis_c(name = "Normalized\nScore") +
    labs(
      title = "Normalized Metric Scores (0-1 Scale)",
      x = "Metric",
      y = "Method"
    ) +
    UtilityTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(paths$figures_dir, "02_heatmap_normalized_metrics.png"), p_heatmap, width = 10, height = 6, bg = "white")

  # --- 4. Radar Charts ---
  log_message("Generating radar charts...")

  radar_data <- league_table %>%
    select(method, starts_with("norm_")) %>%
    rename_with(~gsub("norm_", "", .), starts_with("norm_")) %>%
    tibble::column_to_rownames("method")

  # Add max and min rows for radarchart function
  radar_data <- rbind(rep(1, ncol(radar_data)), rep(0, ncol(radar_data)), radar_data)

  # Generate one chart per method
  for (i in 1:length(config$methods$run)) {
    method_name <- config$methods$run[i]

    png(file.path(paths$figures_dir, paste0("03_radarchart_", method_name, ".png")), width = 800, height = 800, res = 150, bg = "white")
    radarchart(
      df = radar_data[c(1, 2, i + 2), ],
      title = paste("Performance Profile:", method_name),
      pcol = brewer.pal(3, "Set1")[3],
      pfcol = scales::alpha(brewer.pal(3, "Set1")[3], 0.4),
      plwd = 2,
      cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(0, 1, 0.25), cglwd = 0.8,
      vlcex = 0.8
    )
    dev.off()
  }

  # --- 5. UMAP Panels ---
  log_message("Generating UMAP panels...")

  for (method_name in config$methods$run) {
    obj_path <- file.path(analyzed_dir, paste0(method_name, "_analyzed.rds"))
    if (!file.exists(obj_path)) {
      log_message("Skipping UMAP for ", method_name, " (file not found).")
      next
    }
    obj <- readRDS(obj_path)

    umap_reduction_name <- paste0("umap_", method_name)
    batch_var <- config$methods$harmony$group_by_vars[1]
    labels_key <- config$methods$scanvi$labels_key

    p_umap_batch <- DimPlot(obj, reduction = umap_reduction_name, group.by = batch_var, pt.size = 0.1) +
      labs(title = paste(method_name, "- Batch")) + UtilityTheme()
    p_umap_label <- DimPlot(obj, reduction = umap_reduction_name, group.by = labels_key, pt.size = 0.1, label = TRUE, repel = TRUE) +
      labs(title = paste(method_name, "- Cell Type")) + UtilityTheme() + NoLegend()

    p_combined_umap <- p_umap_batch + p_umap_label

    ggsave(file.path(paths$figures_dir, paste0("04_umap_", method_name, ".png")), p_combined_umap, width = 14, height = 6, bg = "white")
  }

  # --- 6. Summary League Table Plot ---
  log_message("Generating summary table plot...")

  table_data <- league_table %>%
    select(Rank, method, GlobalScore, Score_Batch, Score_Label, nmi, ari, ilisi, asw_batch) %>%
    mutate(across(where(is.numeric), ~round(., 3)))

  gt_table <- gt(table_data) %>%
    tab_header(title = "Integration Method Benchmark Rankings") %>%
    cols_label(
      GlobalScore = "Global Score",
      Score_Batch = "Batch Score",
      Score_Label = "Label Score",
      nmi = "NMI",
      ari = "ARI",
      ilisi = "iLISI",
      asw_batch = "1 - ASW"
    ) %>%
    fmt_number(columns = c(GlobalScore, Score_Batch, Score_Label), decimals = 3) %>%
    data_color(
      columns = c(GlobalScore, Score_Batch, Score_Label),
      colors = scales::col_numeric(
        palette = c("white", "cornflowerblue"),
        domain = c(0, 1)
      )
    ) %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = c(Rank, method)))

  gtsave(gt_table, file.path(paths$figures_dir, "05_summary_table.png"))

  log_message("Plot generation finished.")
}

# --- Script Execution ---
if (sys.nframe() == 0) {
  log_message("Running 44_plots.R as a standalone script.")
  # Example:
  # generate_plots("config.yaml")
  log_message("Script 44_plots.R finished.")
}
