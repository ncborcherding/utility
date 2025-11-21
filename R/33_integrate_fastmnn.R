suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(batchelor)
  library(BiocParallel)
  library(BiocNeighbors)      
})

integrate_fastmnn <- function(preprocessed_obj_path, 
                              config_path = "config.yaml") {
  log_message("--- Starting fastMNN Integration ---")
  config <- yaml::read_yaml(config_path)
  obj <- readRDS(preprocessed_obj_path)
  
  integration_timer <- start_timer()
  
  pcs <- obj@reductions$pca@cell.embeddings[, 1:30]
  idx_list <- split(seq_len(nrow(pcs)), obj$orig.ident)
  pcs_list <- lapply(idx_list, function(i) pcs[i, , drop = FALSE])
  rm(pcs)
  
  chunk_size   <- 50         
  k_neighbors  <- config$methods$fastmnn$k
  bnparam      <- AnnoyParam()
  bpparam      <- SerialParam() 
  
  # ---- 1) Make contiguous chunks of batches ----
  n_batches <- length(pcs_list)
  chunks <- split(seq_len(n_batches), ceiling(seq_len(n_batches) / chunk_size))
  
  # ---- 2) Run MNN inside each chunk (Stage 1) ----
  stage1 <- lapply(chunks, function(idx) {
    suppressWarnings(
      reducedMNN(
        pcs_list[idx],
        k           = k_neighbors,
        auto.merge  = TRUE,          
        BNPARAM     = bnparam,
        BPPARAM     = bpparam
      )
    )
  })
  rm(pcs_list)
  
  # ---- 3) Extract corrected embeddings from each chunk ----
  stage1_corr <- lapply(stage1, function(sce) {
    sce$corrected
  })
  rm(stage1)
  
  # ---- 4) Final merge across the chunk-level corrections (Stage 2) ----
  out <- batchelor::reducedMNN(
    stage1_corr,                
    k               = config$k %||% 20,
    prop.k          = config$prop_k %||% NULL,
    ndist           = config$ndist %||% 3,
    min.batch.skip  = config$min_batch_skip %||% 0,
    BNPARAM         = bnparam,
    BPPARAM         = bpparam
  )
  
  fastmnn_embed <- out$corrected[colnames(obj),]
  obj[["fastmnn"]] <- CreateDimReducObject(
    embeddings = fastmnn_embed,
    key = "fastmnn_",
    assay = DefaultAssay(obj)
  )
  
  stop_timer(integration_timer, "fastMNN Integration")
  
  # Join Layers
  obj <- JoinLayers(obj)
  
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "fastmnn.rds")
  safe_save_rds(obj, output_path)
  log_message("Seurat object with fastMNN integration saved.")
  return(output_path)
}
