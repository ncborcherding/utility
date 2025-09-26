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
})

run_cluster_param_grid <- function(
    best_method,
    seurat_rds_path,
    out_dir = "results/07_clustering",
    dims = 1:30,
    graph_name = "KNN",
    resolutions = seq(0.2, 2.0, by = 0.2),
    k_params = c(20, 30, 40, 50),
    weights = NULL,
    random_seed = 42L,
    n_workers = 2
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_message("Loading Seurat object:", seurat_rds_path)
  obj <- readRDS(seurat_rds_path)
  
  grid <- build_param_grid(resolutions, k_params, weights)
  log_message("Total grid size:", nrow(grid))
  
  # Parallel plan
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  
  results <- future_lapply(seq_len(nrow(grid)), function(i) {
    cfg <- grid[i, ]
    obj2 <- obj[[best_method]]
    
    KNN.graph <- bluster::makeKNNGraph(obj[[best_method]]@cell.embeddings,
                                       k = cfg$k_param)
    Cl <- leidenAlg::leiden.community(KNN.graph, 
                                      resolution = cfg$resolution)
    cl_col <- paste0("leiden_", cfg$cfg_id)
    
    obj2[[cl_col]] <- Cl$membership
    obj2[[cl_col]] <- Idents(obj2)
    
    # Metrics (embedding = PCA if available)
    emb <- tryCatch(Embeddings(obj2, reduction = best_method)[, dims, drop = FALSE],
                    error = function(e) NULL)
    sil <- if (!is.null(emb)) compute_silhouette(emb, obj2[[cl_col, drop = TRUE]]) else NA_real_
    mod <- compute_modularity(KNN.graph, obj2, cl_col)
    conn <- compute_graph_connectivity(KNN.graph, obj2, cl_col)
    
    tibble::tibble(
      cfg_id = cfg$cfg_id,
      resolution = cfg$resolution,
      k_param = cfg$k_param,
      prune_snn = cfg$prune_snn,
      silhouette = sil,
      modularity = mod,
      connectivity = conn,
      n_clusters = nlevels(obj2[[cl_col, drop = TRUE]]),
      n_singletons = sum(table(obj2[[cl_col, drop = TRUE]]) == 1)
    )
  }) %>% bind_rows()
  
  readr::write_csv(results, file.path(out_dir, "grid_metrics.csv"))
  saveRDS(results, file.path(out_dir, "grid_metrics.rds"))
  log_message("Saved metrics to:", file.path(out_dir, "grid_metrics.csv"))
  
  invisible(results)
}
