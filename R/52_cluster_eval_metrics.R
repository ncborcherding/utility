# R/52_cluster_eval_metrics.R
# Combine grid metrics and stability; compute a composite score

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(scales)
})

# Robust min-max rescale (winsorize extremes)
.rescale01 <- function(x) {
  if (all(is.na(x))) return(x)
  q <- quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)
  x <- pmin(pmax(x, q[1]), q[2])
  scales::rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE))
}

combine_and_score <- function(grid_metrics_path = "results/07_clustering/grid_metrics.rds",
                              stability_path     = "results/07_clustering/grid_stability.rds",
                              out_csv            = "results/07_clustering/grid_scores.csv",
                              w_sil = 0.25, 
                              w_mod = 0.25, 
                              w_conn = 0.25, 
                              w_stab = 0.25,
                              singleton_penalty = 0.05) {
  grid <- readRDS(grid_metrics_path)
  stab <- readRDS(stability_path)
  all <- dplyr::left_join(grid, 
                          stab, 
                          by = c("cfg_id", "resolution", "k_param"))
  
  all <- all %>% mutate(
    z_sil  = .rescale01(silhouette),
    z_mod  = .rescale01(modularity),
    z_conn = .rescale01(connectivity),
    z_stab = .rescale01(stability_jaccard)
  ) %>% mutate(
    composite = w_sil * z_sil + w_mod * z_mod + w_conn * z_conn + w_stab * z_stab -
      singleton_penalty * (n_singletons > 0)
  ) %>% arrange(desc(composite))
  
  readr::write_csv(all, out_csv)
  saveRDS(all, sub("\\.csv$", ".rds", out_csv))
}
