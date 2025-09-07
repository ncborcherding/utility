# R/51_cluster_stability.R


# Subsampled fits and Jaccard stability vs full
n <- ncol(obj)
idx_all <- colnames(obj)
jaccs <- replicate(n_repeats, {
  idx <- sort(sample.int(n, size = floor(subsample_frac * n), replace = FALSE))
  cells <- idx_all[idx]
  obj_sub <- subset(obj, cells = cells)
  obj_sub <- safe_FindNeighbors(obj_sub, dims = dims, k_param = cfg$k_param,
                                graph_name = graph_name, prune.SNN = cfg$prune_snn)
  obj_sub <- safe_FindClusters(obj_sub, graph_name = graph_name, resolution = cfg$resolution,
                               algorithm = 4L, random_seed = seed)
  cl_sub <- as.character(Idents(obj_sub))
  
  
  # Compare on intersection of cells
  common <- intersect(names(cl_full), names(cl_sub))
  if (length(common) < 10) return(NA_real_)
  jaccard_stability(cl_full[common], cl_sub[common])
})


tibble::tibble(
  cfg_id = cfg$cfg_id,
  resolution = cfg$resolution,
  k_param = cfg$k_param,
  prune_snn = cfg$prune_snn,
  stability_jaccard = mean(jaccs, na.rm = TRUE)
)
}


run_cluster_stability <- function(
    seurat_rds_path,
    grid_rds_path = "results/05_clustering/grid_metrics.rds",
    out_dir = "results/05_clustering",
    dims = 1:30,
    graph_name = "SNN",
    subsample_frac = 0.8,
    n_repeats = 25,
    n_workers = 4,
    seed = 1L
) {
  ensure_leiden_available()
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  obj <- readRDS(seurat_rds_path)
  grid <- readRDS(grid_rds_path) %>% distinct(cfg_id, resolution, k_param, prune_snn)
  
  
  oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  
  
  # Convert grid rows to list of cfgs
  cfgs <- split(grid, seq_len(nrow(grid)))
  
  
  stab <- future_lapply(cfgs, function(cfg) {
    stability_for_config(obj, dims, graph_name, cfg[[1]],
                         subsample_frac = subsample_frac,
                         n_repeats = n_repeats, seed = seed)
  }) %>% bind_rows()
  
  
  readr::write_csv(stab, file.path(out_dir, "grid_stability.csv"))
  saveRDS(stab, file.path(out_dir, "grid_stability.rds"))
  log_message("Saved stability metrics to:", file.path(out_dir, "grid_stability.csv"))
  
  
  invisible(stab)
}


if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript R/51_cluster_stability.R <seurat_rds>")
  run_cluster_stability(args[[1]])
}