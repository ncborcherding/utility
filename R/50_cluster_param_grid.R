# R/50_cluster_param_grid.R
# Sweep Leiden parameters on an integrated Seurat object
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(bluster)
  library(leidenAlg)
  library(future.apply)
})

# ---- helpers you already have ----
# build_param_grid(), compute_silhouette(), compute_modularity(), compute_graph_connectivity()
# log_message()

.run_one_cfg <- function(cfg, seurat_rds_path, best_method, dims) {
  obj <- readRDS(seurat_rds_path)
  
  # use embedding matrix only; no need to carry full Seurat around for metrics
  emb <- tryCatch(
    Embeddings(obj, reduction = best_method)[, dims, drop = FALSE],
    error = function(e) NULL
  )
  if (is.null(emb)) {
    stop(sprintf("No embeddings found for reduction '%s'", best_method))
  }
  
  # build graph on the embedding (equivalent to obj[[best_method]]@cell.embeddings)
  g <- bluster::makeKNNGraph(emb, k = cfg$k_param)
  
  cl <- leidenAlg::leiden.community(g, resolution = cfg$resolution)
  membership <- cl$membership
  
  sil  <- compute_silhouette(emb, membership)
  mod  <- compute_modularity(g, membership)         
  conn <- compute_graph_connectivity(g, membership)  
  
  tibble::tibble(
    cfg_id       = cfg$cfg_id,
    resolution   = cfg$resolution,
    k_param      = cfg$k_param,
    prune_snn    = cfg$prune_snn %||% NA,
    silhouette   = sil,
    modularity   = mod,
    connectivity = conn,
    n_clusters   = length(unique(membership)),
    n_singletons = sum(tabulate(membership)[tabulate(membership) == 1])
  )
}

run_cluster_param_grid <- function(
    best_method,
    seurat_rds_path,
    out_dir      = "results/07_clustering",
    dims         = 1:30,
    graph_name   = "KNN",
    resolutions  = seq(0.2, 2.0, by = 0.2),
    k_params     = c(20, 30, 40, 50),
    weights      = NULL,
    random_seed  = 42L,
    n_workers    = 2
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_message("Loading Seurat object (once, to build grid):", seurat_rds_path)
  obj <- readRDS(seurat_rds_path)
  
  grid <- build_param_grid(resolutions, k_params, weights)
  log_message("Total grid size:", nrow(grid))
  
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  
  set.seed(random_seed)
  
  results <- future.apply::future_lapply(
    seq_len(nrow(grid)),
    FUN = function(i) .run_one_cfg(grid[i, ], seurat_rds_path, best_method, dims),
    future.seed = TRUE,
    future.globals = c(".run_one_cfg", "grid", "seurat_rds_path", "best_method", "dims")
  ) %>% dplyr::bind_rows()
  
  readr::write_csv(results, file.path(out_dir, "grid_metrics.csv"))
  saveRDS(results, file.path(out_dir, "grid_metrics.rds"))
  log_message("Saved metrics to:", file.path(out_dir, "grid_metrics.csv"))
  
  invisible(results)
}