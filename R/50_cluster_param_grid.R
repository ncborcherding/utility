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
  library(SeuratObject)
  library(future.apply)
})

run_cluster_param_grid <- function(best_method,
                                   seurat_rds_path,
                                   out_dir      = "results/07_clustering",
                                   dims         = 1:30,
                                   graph_name   = "KNN",
                                   resolutions  = seq(0.2, 2.0, by = 0.2),
                                   k_params     = c(10, 20, 30, 40, 50),
                                   weights      = NULL,
                                   random_seed  = 42L,
                                   n_workers    = 2) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_message("Loading Seurat object (once, to build grid):", seurat_rds_path)
  obj <- readRDS(seurat_rds_path)
  emb <- Embeddings(obj, reduction = best_method)[, dims, drop = FALSE]
  rm(obj); gc()
  
  grid <- build_param_grid(resolutions, k_params, weights)
  log_message("Total grid size:", nrow(grid))
  
  # split work into N chunks (one per worker)
  n <- NROW(grid)  # works for data.frames, matrices, and vectors
  idx_chunks <- split(
    seq_len(n),
    rep_len(seq_len(n_workers), n)
  )
  
  
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  
  set.seed(random_seed)
  
  chunk_results <- future.apply::future_lapply(
    idx_chunks,
    function(idxs) {
      out <- vector("list", length(idxs))
      for (j in seq_along(idxs)) {
        i <- idxs[j]
        cfg <- grid[i, ]
        g <- bluster::makeKNNGraph(emb, k = cfg$k_param)
        cl <- leidenAlg::leiden.community(g, resolution = cfg$resolution)
        membership <- cl$membership
        
        out[[j]] <- tibble::tibble(
          cfg_id       = cfg$cfg_id,
          resolution   = cfg$resolution,
          k_param      = cfg$k_param,
          silhouette   = compute_silhouette(emb, membership),
          modularity   = compute_modularity(g, membership),
          connectivity = compute_graph_connectivity(g, membership),
          n_clusters   = length(unique(membership)),
          n_singletons = sum(tabulate(membership) == 1)
        )
      }
      dplyr::bind_rows(out)
    },
    future.globals = c("grid","emb","compute_silhouette",
                       "compute_modularity","compute_graph_connectivity")
  )
  
  results <- dplyr::bind_rows(chunk_results)
  
  readr::write_csv(results, file.path(out_dir, "grid_metrics.csv"))
  saveRDS(results, file.path(out_dir, "grid_metrics.rds"))
  log_message("Saved metrics to:", file.path(out_dir, "grid_metrics.csv"))
  
  invisible(results)
}
