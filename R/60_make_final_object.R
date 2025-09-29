# R/60_make_final_object.R

suppressPackageStartupMessages({
  library(Seurat)
  library(bluster)
  library(leidenAlg)
})

make_final_object <- function(cfg, 
                              best_method, 
                              k, 
                              resolution) {
  
  integration_func <- get(paste0("integrate_", best_method))
  integrated_obj_path <- integration_func(preprocessed_obj_path, config_path)
  
  
  #TODO May need to refine this as > 1 million cells is computationally massive
  g <- bluster::makeKNNGraph(emb, k = k)
  cl <- leidenAlg::leiden.community(g, resolution = resolution)
  obj$leiden_clusters <- cl$membership
  Idents(obj) <- "leiden_clusters"
}