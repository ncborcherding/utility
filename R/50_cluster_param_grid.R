# R/50_cluster_param_grid.R
log_message("Loading Seurat object:", seurat_rds_path)
obj <- readRDS(seurat_rds_path)


grid <- build_param_grid(resolutions, k_params, prune_snn, weights)
log_message("Total grid size:", nrow(grid))


# Parallel plan
oplan <- future::plan()
on.exit(future::plan(oplan), add = TRUE)
future::plan(future::multisession, workers = n_workers)


results <- future_lapply(seq_len(nrow(grid)), function(i) {
  cfg <- grid[i, ]
  SeuratObject::Idents(obj) <- NULL
  obj2 <- obj
  # 1) Neighbors
  obj2 <- safe_FindNeighbors(
    obj2, dims = dims, k_param = cfg$k_param, graph_name = graph_name,
    prune.SNN = cfg$prune_snn
  )
  # 2) Leiden
  obj2 <- safe_FindClusters(
    obj2, graph_name = graph_name, resolution = cfg$resolution,
    algorithm = 4L, random_seed = random_seed
  )
  cl_col <- paste0("leiden_", cfg$cfg_id)
  obj2[[cl_col]] <- Idents(obj2)
  
  
  # Metrics (embedding = PCA if available)
  emb <- tryCatch(Embeddings(obj2, reduction = "pca")[, dims, drop = FALSE],
                  error = function(e) NULL)
  sil <- if (!is.null(emb)) compute_silhouette(emb, obj2[[cl_col, drop = TRUE]]) else NA_real_
  mod <- compute_modularity(obj2, cl_col, graph_name)
  conn <- compute_graph_connectivity(obj2, cl_col, graph_name)
  
  
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


if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript R/50_cluster_param_grid.R <seurat_rds>")
  run_cluster_param_grid(args[[1]])
}