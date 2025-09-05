suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(yaml)
  library(batchelor)
  library(BiocParallel)
  library(BiocNeighbors)      # <- add this
  suppressWarnings(suppressMessages(requireNamespace("RhpcBLASctl", quietly = TRUE)))
})

cap_threads <- function(blas_threads = 1) {
  # Pin BLAS/OpenMP threads to avoid oversubscription
  Sys.setenv(
    OMP_NUM_THREADS = blas_threads,
    OPENBLAS_NUM_THREADS = blas_threads,
    MKL_NUM_THREADS = blas_threads,
    VECLIB_MAXIMUM_THREADS = blas_threads,
    NUMEXPR_NUM_THREADS = blas_threads
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(blas_threads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(blas_threads), silent = TRUE)
  }
}

mk_bpparam <- function(workers = 8) {
  workers <- max(1L, as.integer(workers))
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = TRUE)
  } else {
    BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
  }
}

integrate_fastmnn <- function(preprocessed_obj_path, config_path = "config.yaml") {
  log_message("--- Starting fastMNN (on PCA) Integration ---")
  config <- yaml::read_yaml(config_path)
  method_config <- config$methods$fastmnn
  obj <- readRDS(preprocessed_obj_path)
  set.seed(config$seed)
  
  # ----- Thread & parallel caps -----
  # Read desired workers from config if present; default to 8
  workers <- method_config$workers %||% 8
  cap_threads(blas_threads = 1)
  BPPARAM <- mk_bpparam(workers = workers)
  log_message("Using ", class(BPPARAM)[1], " with ", workers, " worker(s).")
  
  # ----- PCA checks/truncation -----
  if (!"pca" %in% Reductions(obj)) stop("PCA reduction not found; run PCA first.")
  pca_embed <- Embeddings(obj, "pca")
  n_use <- min(ncol(pca_embed), method_config$n_dims %||% 40L)  # <= modest cap
  if (n_use < 2) stop("Need at least 2 PCs; got n_dims = ", n_use)
  pca_embed <- pca_embed[, seq_len(n_use), drop = FALSE]
  
  # ----- Batch factor -----
  meta <- obj@meta.data
  batch_var <- method_config$batch_var
  if (!batch_var %in% colnames(meta)) {
    stop("fastMNN 'batch_var' (", batch_var, ") not found in metadata.")
  }
  batch <- factor(meta[[batch_var]])
  log_message("Running reducedMNN on ", n_use, " PCs across ", nlevels(batch), " batch(es).")
  
  # ----- Fast approximate NN index -----
  # HNSW is very fast & accurate for high-dimensional PCs
  bnparam <- BiocNeighbors::HnswParam(nlinks = method_config$hnsw_M %||% 16,
                                      ef.search = method_config$hnsw_ef_search %||% 50)
  
  # ----- Run reducedMNN with safe merge & capped parallelism -----
  integration_timer <- start_timer()
  mnn_out <- batchelor::reducedMNN(
    pca_embed,
    batch           = batch,
    k               = method_config$k %||% 20,
    prop.k          = method_config$prop_k %||% NULL,
    ndist           = method_config$ndist %||% 3,
    merge.order     = NULL,          # let auto.merge design a safe plan
    auto.merge      = TRUE,
    min.batch.skip  = method_config$min_batch_skip %||% 0,
    BNPARAM         = bnparam,       # <- approximate neighbors
    BPPARAM         = BPPARAM        # <- capped workers
  )
  
  fastmnn_embed <- mnn_out$corrected
  rownames(fastmnn_embed) <- rownames(pca_embed)
  obj[["fastmnn"]] <- CreateDimReducObject(
    embeddings = fastmnn_embed,
    key = "fastmnn_",
    assay = DefaultAssay(obj)
  )
  
  stop_timer(integration_timer, "fastMNN (on PCA) Integration")
  
  output_dir <- file.path(config$paths$results_dir, "04_integrated_data")
  output_path <- file.path(output_dir, "fastmnn.rds")
  safe_save_rds(obj, output_path)
  log_message("Seurat object with fastMNN (on PCA) integration saved.")
  return(output_path)
}
