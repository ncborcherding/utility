# R/00_clust_utils.R
# Helper utilities for parameter grids, Leiden checks, and metrics

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(future.apply)
  library(igraph)
})

# --- Logging ---
log_message <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(...)))
}

# --- Build parameter grid ---
build_param_grid <- function(resolutions, k_params, weights = NULL) {
  expand.grid(
    resolution = resolutions,
    k_param    = k_params,
    weights    = if (is.null(weights)) NA_character_ else weights,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    mutate(cfg_id = sprintf("res%.2f_k%d%s",
                            resolution, k_param,
                            ifelse(is.na(weights), "", paste0("_w", weights))))
}

# --- Silhouette on an embedding ---
# Uses PCA by default. Provide a matrix 'emb' (cells x dims) and
# factor 'cl' of cluster labels.
compute_silhouette <- function(emb, cl) {
  suppressWarnings({
    d <- stats::dist(emb)
    sil <- cluster::silhouette(as.integer(cl), d)
  })
  mean(sil[, "sil_width"])
}

# --- Graph modularity for the cluster partition ---
compute_modularity <- function(g, membership) {
  if (is.null(g) || is.null(membership)) return(NA_real_)
  if (igraph::vcount(g) != length(membership)) return(NA_real_)
  cl_int <- as.integer(factor(membership))  # robust to non-contiguous labels
  igraph::modularity(g, membership = cl_int)
}

# --- Graph connectivity: average fraction of nodes in LCC within clusters ---
# (Higher is more connected; 1.0 means each cluster is a single connected piece.)
compute_graph_connectivity <- function(g, membership) {
  if (is.null(g) || is.null(membership)) return(NA_real_)
  if (igraph::vcount(g) != length(membership)) return(NA_real_)
  
  cl <- factor(membership)
  # Split node indices by cluster (fast; avoids repeated which())
  idx_list <- split(seq_along(cl), cl)
  
  fracs <- vapply(idx_list, function(nodes) {
    n <- length(nodes)
    if (n <= 1) return(1)  # singleton cluster treated as fully connected
    sg <- igraph::induced_subgraph(g, vids = nodes)
    comps <- igraph::components(sg)$csize
    max(comps) / n
  }, numeric(1))
  
  mean(fracs)
}



# --- Jaccard stability between two labelings ---
jaccard_stability <- function(cl1, cl2) {
  stopifnot(length(cl1) == length(cl2))
  cl1 <- factor(cl1); cl2 <- factor(cl2)
  tab <- table(cl1, cl2)
  jacc <- outer(levels(cl1), levels(cl2), Vectorize(function(a, b) {
    inter <- sum(cl1 == a & cl2 == b)
    u <- sum(cl1 == a) + sum(cl2 == b) - inter
    if (u == 0) return(1)
    inter / u
  }))
  if(dim(jacc)[1] > dim(jacc)[2]) {
    jacc <- t(jacc)
  }
  # maximize sum via Hungarian (solve assignment on cost = 1 - jacc)
  assignment <- clue::solve_LSAP(1 - jacc)
  mean(jacc[cbind(seq_along(assignment), assignment)])
}
