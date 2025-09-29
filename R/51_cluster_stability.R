# R/51_cluster_stability.R
# Subsampling stability assessment across parameter grid

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(future.apply)
  library(Seurat)
  library(bluster)
  library(leidenAlg)
  library(clue)   
  library(igraph)
})

# --- One clustering pass using the SAME method as above (bluster + leidenAlg) ---
cluster_once_leiden <- function(emb, 
                                k, 
                                resolution, 
                                seed = NULL, 
                                directed = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  g  <- bluster::makeKNNGraph(emb, k = k, directed = directed)
  cl <- leidenAlg::leiden.community(g, resolution = resolution)
  memb <- cl$membership
  names(memb) <- rownames(emb)      # keep cell names aligned
  list(membership = memb, graph = g)
}

# --- Stability for a single config vs subsamples, using the same clustering ---
stability_for_config_emb <- function(emb,                 
                                     cfg,                 
                                     subsample_frac = 0.8,
                                     n_repeats      = 25,
                                     seed           = 42L) {
  stopifnot(is.matrix(emb) || inherits(emb, "Matrix"))
  stopifnot(!is.null(rownames(emb)))   # must have cell names
  
  n <- nrow(emb)
  if (n < 2L) {
    return(tibble::tibble(
      cfg_id = cfg$cfg_id,
      resolution = cfg$resolution,
      k_param = cfg$k_param,
      stability_jaccard = NA_real_
    ))
  }
  
  # Full fit
  full <- cluster_once_leiden(emb, 
                              k = cfg$k_param, 
                              resolution = cfg$resolution, 
                              seed = seed)
  cl_full <- full$membership
  
  # Repeats on subsamples
  set.seed(seed)
  sizesub <- max(2L, floor(subsample_frac * n))
  jaccs <- numeric(n_repeats)
  
  for (r in seq_len(n_repeats)) {
    
    idx  <- sort(sample.int(n, size = sizesub, replace = FALSE))
    sub_emb <- emb[idx, , drop = FALSE]
    
    sub  <- cluster_once_leiden(sub_emb, k = cfg$k_param, resolution = cfg$resolution)
    cl_sub <- sub$membership
    
    common <- names(cl_sub)                    # subsample cells
    if (length(common) < 10L) {
      jaccs[r] <- NA_real_
    } else {
      jaccs[r] <- jaccard_stability(cl_full[common], cl_sub[common])
    }
  }
  
  tibble::tibble(
    cfg_id            = cfg$cfg_id,
    resolution        = cfg$resolution,
    k_param           = cfg$k_param,
    stability_jaccard = mean(jaccs, na.rm = TRUE)
  )
}

run_cluster_stability <- function(seurat_rds_path,
                                  best_method,
                                  grid_rds_path = "results/07_clustering/grid_metrics.rds",
                                  out_dir = "results/07_clustering",
                                  dims = 1:30,
                                  subsample_frac = 0.8,
                                  n_repeats = 25,
                                  n_workers = 2,
                                  seed = 42L) {

  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- Load once, keep only the embedding ---
  obj <- readRDS(seurat_rds_path)
  emb_all <- Embeddings(obj, reduction = best_method)
  if (max(dims) > ncol(emb_all)) stop("Requested dims exceed available dimensions.")
  emb <- emb_all[, dims, drop = FALSE]
  if (is.null(rownames(emb))) rownames(emb) <- colnames(obj)
  rm(obj, emb_all); gc()
  
  # --- Read the grid (unique cfgs) ---
  grid <- readRDS(grid_rds_path) %>% dplyr::distinct(cfg_id, resolution, k_param)
  if (nrow(grid) == 0L) stop("No grid rows found in grid_rds_path.")
  
  # Convert rows to plain named lists so $access works cleanly
  cfgs <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
  
  # --- Parallel plan ---
  oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  
  # Chunk configs so each worker receives 'emb' only once
  idx_chunks <- split(seq_along(cfgs), rep_len(seq_len(n_workers), length(cfgs)))
  
  chunk_results <- future.apply::future_lapply(
    idx_chunks,
    function(idxs) {
      out <- vector("list", length(idxs))
      for (j in seq_along(idxs)) {
        cfg <- cfgs[[ idxs[j] ]]  # cfg is a named list with cfg_id, resolution, k_param
        out[[j]] <- stability_for_config_emb(
          emb,
          cfg,
          subsample_frac = subsample_frac,
          n_repeats      = n_repeats,
          seed           = seed
        )
      }
      dplyr::bind_rows(out)
    },
    future.seed    = TRUE,
    future.globals = c("cfgs","emb",
                       "stability_for_config_emb",
                       "cluster_once_leiden",
                       "jaccard_stability",
                       "subsample_frac","n_repeats","seed")
  )
  
  stab <- dplyr::bind_rows(chunk_results)
  
  readr::write_csv(stab, file.path(out_dir, "grid_stability.csv"))
  saveRDS(stab, file.path(out_dir, "grid_stability.rds"))
  if (exists("log_message")) log_message("Saved stability metrics to:", file.path(out_dir, "grid_stability.csv"))
  
  invisible(stab)
}