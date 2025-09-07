# R/00_clust_utils.R
obj <- Seurat::FindClusters(
  object = obj, graph.name = graph_name,
  resolution = resolution, algorithm = algorithm, ...
)
obj
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
compute_modularity <- function(obj, cluster_col, graph_name = "SNN") {
  g <- obj@graphs[[graph_name]]
  if (is.null(g)) return(NA_real_)
  if (!igraph::is.igraph(g)) g <- SeuratObject::as.graph(g)
  cl <- factor(obj[[cluster_col, drop = TRUE]])
  igraph::modularity(g, membership = as.integer(cl))
}


# --- Graph connectivity: average fraction of nodes in LCC within clusters ---
# (Higher is more connected; 1.0 means each cluster is a single connected piece.)
compute_graph_connectivity <- function(obj, cluster_col, graph_name = "SNN") {
  g <- obj@graphs[[graph_name]]
  if (is.null(g)) return(NA_real_)
  if (!igraph::is.igraph(g)) g <- SeuratObject::as.graph(g)
  cl <- factor(obj[[cluster_col, drop = TRUE]])
  fracs <- vapply(levels(cl), function(cc) {
    nodes <- which(cl == cc)
    if (length(nodes) <= 1) return(1)
    sg <- igraph::induced_subgraph(g, vids = nodes)
    comps <- igraph::components(sg)$csize
    max(comps) / sum(comps)
  }, numeric(1))
  mean(fracs)
}


# --- Jaccard stability between two labelings ---
# Computes average Jaccard index between matched clusters using Hungarian algorithm
jaccard_stability <- function(cl1, cl2) {
  cl1 <- factor(cl1); cl2 <- factor(cl2)
  # contingency table
  tab <- table(cl1, cl2)
  # convert to Jaccard matrix
  jacc <- outer(levels(cl1), levels(cl2), Vectorize(function(a, b) {
    inter <- sum(cl1 == a & cl2 == b)
    u <- sum(cl1 == a) + sum(cl2 == b) - inter
    if (u == 0) return(1)
    inter / u
  }))
  # maximize sum via Hungarian (solve assignment on cost = 1 - jacc)
  assignment <- clue::solve_LSAP(1 - jacc)
  mean(jacc[cbind(seq_along(assignment), assignment)])
}