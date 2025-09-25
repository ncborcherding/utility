# R/53_choose_clusters.R
# Apply the selected configuration, save final clusters into Seurat object


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})


source("R/00_clust_utils.R")


apply_best_clusters <- function(
    seurat_rds_in,
    scores_csv = "results/05_clustering/grid_scores.csv",
    seurat_rds_out = "results/05_clustering/final_clustered_object.rds",
    dims = 1:30,
    graph_name = "SNN",
    top_n = 1L,
    seed = 1L
) {
  ensure_leiden_available()
  obj <- readRDS(seurat_rds_in)
  scores <- readr::read_csv(scores_csv, show_col_types = FALSE)
  best <- dplyr::slice_max(scores, order_by = composite, n = top_n)
  
  
  # Use the single best
  cfg <- best[1, ]
  log_message("Applying best cfg:", cfg$cfg_id)
  
  
  obj <- safe_FindNeighbors(obj, dims = dims, k_param = cfg$k_param, graph_name = graph_name,
                            prune.SNN = cfg$prune_snn)
  obj <- safe_FindClusters(obj, graph_name = graph_name, resolution = cfg$resolution,
                           algorithm = 4L, random_seed = seed)
  
  
  final_col <- paste0("leiden_final_", cfg$cfg_id)
  obj[[final_col]] <- Idents(obj)
  
  
  saveRDS(obj, seurat_rds_out)
  # Also write a light-weight per-cell table for downstream viz
  cell_df <- data.frame(cell = colnames(obj), cluster = obj[[final_col, drop = TRUE]])
  readr::write_csv(cell_df, file.path(dirname(seurat_rds_out), "final_clusters.csv"))
  
  
  invisible(list(object_path = seurat_rds_out, final_label = final_col))
}


if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript R/53_choose_clusters.R <seurat_rds_in>")
  apply_best_clusters(args[[1]])
}