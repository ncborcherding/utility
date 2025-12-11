# R/60_make_final_object.R

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(yaml)
  library(future)
  library(future.apply)
  library(scran)
  library(scRepertoire)
})

source("R/00_utils.R")

#' Generate the Full BPCells Directory
generate_full_bpcells <- function(config) {
  log_message("--- Generating FULL BPCells Dataset (Skipping Subset) ---")
  
  paths <- config$paths
  filter_config <- config$filtering
  
  # Output directory for FULL BPCells matrices
  bp_dir_full <- file.path(paths$results_dir, "01_bpcells_data_FULL")
  
  # If it exists, we assume it's done
  if (dir.exists(bp_dir_full)) {
    log_message("Full BPCells directory already exists. Using existing data.")
    return(bp_dir_full)
  }
  dir.create(bp_dir_full, recursive = TRUE)
  
  obj_dir <- paths$seurat_objects
  rds_files <- list.files(obj_dir, pattern = "\\.rds$", full.names = TRUE)
  
  for (file in rds_files) {
    sample_name <- tools::file_path_sans_ext(basename(file))
    
    obj <- readRDS(file)
    
    # Apply Filtering (Same as subset run)
    if (filter_config$run) {
      meta <- obj@meta.data
      cells_to_keep <- rep(TRUE, ncol(obj))
      for (f in filter_config$filters) {
        if (f$column %in% colnames(meta)) {
          col_data <- meta[[f$column]]
          if (f$type == "in") cells_to_keep <- cells_to_keep & (col_data %in% f$values)
          else if (f$type == "equals") cells_to_keep <- cells_to_keep & (col_data == f$values)
          else if (f$type == "not_na") cells_to_keep <- cells_to_keep & !is.na(col_data)
        }
      }
      if (sum(cells_to_keep) == 0) next
      obj <- subset(obj, cells = colnames(obj)[which(cells_to_keep)])
    }
    
    if (ncol(obj) == 0) next
    
    sample_bp_dir <- file.path(bp_dir_full, sample_name)
    
    counts_mat <- obj@assays$RNA@layers$counts
    if (is.null(counts_mat)) counts_mat <- obj@assays$RNA@counts
    
    BPCells::write_matrix_dir(mat = counts_mat, dir = sample_bp_dir)
    
    # Save metadata/features
    saveRDS(obj[[]], file = file.path(bp_dir_full, paste0(sample_name, "_meta.rds")))
    saveRDS(rownames(obj), file = file.path(bp_dir_full, paste0(sample_name, "_features.rds")))
    
    rm(obj, counts_mat); gc()
  }
  
  return(bp_dir_full)
}

make_final_object <- function(cfg, best_method, k, resolution) {
  
  log_message("=== Starting Creation of Final Full Object ===")
  log_message("Selected Best Method: ", best_method)
  log_message("Selected K: ", k, " | Resolution: ", resolution)
  
  # 1. Ensure Full Data Exists
  full_bp_dir <- generate_full_bpcells(cfg)
  
  # 2. Load Full Data
  log_message("Loading Full BPCells matrices...")
  sample_dirs <- list.dirs(full_bp_dir, full.names = FALSE, recursive = FALSE)
  sample_names <- sample_dirs[!grepl("_meta.rds|_features.rds", sample_dirs)]
  
  # Set up parallel plan for reading/processing
  n_workers <- cfg$memory$num_workers %||% 4
  if(n_workers == 0) n_workers <- future::availableCores() - 12
  plan("multisession", workers = n_workers)
  
  # Load objects into a list (lightweight pointers)
  object_list <- lapply(sample_names, function(s_name) {
    counts_bp <- BPCells::open_matrix_dir(file.path(full_bp_dir, s_name))
    meta <- readRDS(file.path(full_bp_dir, paste0(s_name, "_meta.rds")))
    rownames(counts_bp) <- readRDS(file.path(full_bp_dir, paste0(s_name, "_features.rds")))
    colnames(counts_bp) <- rownames(meta)
    Seurat::CreateSeuratObject(counts = counts_bp, meta.data = meta)
  })
  
  # 3. Preprocessing: Batch-aware HVG Calculation (scran style)
  log_message("Calculating HVGs per sample using scran::modelGeneVar (Parallel)...")
  
  # Use future_lapply for speed on 700+ samples
  per_sample_stats <- future_lapply(object_list, function(obj) {
    # Access BPCells counts
    counts_bp <- obj@assays$RNA@layers$counts
    if (is.null(counts_bp)) counts_bp <- obj@assays$RNA@counts
    
    # Calculate library sizes efficiently
    lib_sizes <- as.numeric(BPCells::colSums(counts_bp))
    lib_sizes[lib_sizes == 0] <- 1
    scale_factor <- median(lib_sizes)
    
    # Normalize (lazy op)
    norm_bp <- t(t(counts_bp) / lib_sizes) * scale_factor
    
    norm_mat <- as(norm_bp, "dgCMatrix")
    logcounts <- log1p(norm_mat)
    
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = logcounts))
    diff <- scran::modelGeneVar(sce, assay.type = "logcounts")
    diff@rownames <- rownames(obj)
    return(diff)
  }, future.seed = TRUE)
  
  # Combine variance stats
  log_message("Combining variance statistics...")
  common_genes <- Reduce(intersect, lapply(per_sample_stats, rownames))
  per_sample_stats <- lapply(per_sample_stats, function(stat) stat[common_genes, , drop = FALSE])
  combined_stats <- scran::combineVar(per_sample_stats)
  
  n_features <- cfg$preprocess$n_variable_features
  top_hvg <- scran::getTopHVGs(combined_stats, n = n_features)
  
  # Filter unwanted genes (MT, Ribosomal, VDJ)
  variable_features <- scRepertoire::quietVDJgenes(top_hvg)
  variable_features <- variable_features[!grepl('^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-', variable_features)]
  log_message("Identified ", length(variable_features), " robust variable features.")
  
  # 4. Merge and Scale
  log_message("Merging into single BPCells object...")
  if (length(object_list) > 1) {
    obj.full <- merge(object_list[[1]], y = object_list[2:length(object_list)])
  } else {
    obj.full <- object_list[[1]]
  }
  rm(object_list); gc()
  
  VariableFeatures(obj.full) <- variable_features
  
  log_message("Scaling and Running PCA on Full Dataset...")
  obj.full <- NormalizeData(obj.full, verbose = FALSE)
  obj.full <- ScaleData(obj.full, features = variable_features, verbose = FALSE)
  obj.full <- RunPCA(obj.full, npcs = 50, verbose = FALSE)
  
  # 5. Integration
  log_message("Running Integration: ", best_method)
  
  # Save preprocessed full object if needed by integration scripts
  full_preprocessed_path <- file.path(cfg$paths$results_dir, "03_preprocessed_data", "full_object_preprocessed.rds")
  safe_save_rds(obj.full, full_preprocessed_path)
  
  integration_func <- get(paste0("integrate_", best_method))
  full_integrated_path <- integration_func(full_preprocessed_path, "config.yaml")
  
  # Reload for clustering
  obj.integrated <- readRDS(full_integrated_path)
  
  # 6. UMAP Embedding
  reduction_use <- if(best_method == "harmony") "harmony" else 
    if(best_method %in% c("scvi", "scanvi")) "integrated.dr" else "pca"
  if(!reduction_use %in% names(obj.integrated@reductions)) reduction_use <- "pca"
  
  obj.integrated <- RunUMAP(obj.integrated,
                            reduction = reduction_use,
                            dims = 1:config$post_integration$n_dims_use,
                            verbose = FALSE)

  # 7. Clustering 
  emb <- Embeddings(obj.integrated, reduction = reduction_use)
  cluster_composite_score <- combine_and_score()
  
  log_message("Generating Nearest Neighbors (k = ", cluster_composite_score$k, ")...")
  g <- bluster::makeKNNGraph(emb, k = cluster_composite_score$k)
  
  log_message("Running Leiden Clustering (Resolution ", cluster_composite_score$resolution, ")...")
  cl <- leidenAlg::leiden.community(g, resolution = cluster_composite_score$resolution)
  
  obj.integrated$leiden.cluster <- cl$membership
  Idents(obj.integrated) <- "leiden.cluster"
  
  # 8. Generate Plots
  generate_final_plots(obj.integrated, cfg, output_prefix = "06_")

  # 9. Save Final
  final_out_path <- file.path(cfg$paths$results_dir, "FINAL_integrated_object.rds")
  log_message("Saving Final Object to: ", final_out_path)
  safe_save_rds(obj.integrated, final_out_path)
  
  log_message("=== Final Object Created Successfully ===")
  return(final_out_path)
}

  # =============================================================================
  # HELPER: Generate Summary Plots
  # =============================================================================

#' Generate summary UMAP plots for the final integrated object
generate_final_plots <- function(obj, config, output_prefix = "FINAL") {
  
  log_message("--- Generating Summary Plots ---")
  
  fig_dir <- config$paths$figures_dir
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  batch_var <- config$methods$harmony$group_by_vars[1]
  labels_l1 <- "predicted.celltype.l1"
  labels_l2 <- config$post_integration$labels %||% "predicted.celltype.l2"
  
  available_vars <- colnames(obj@meta.data)
  
  # Use rasterization for large datasets
  use_raster <- ncol(obj) > 50000
  raster_dpi <- c(512, 512)
  pt_size <- if(ncol(obj) > 500000) 0.8 else if(ncol(obj) > 100000) 1.6 else 2
  
  # ---------------------------------------------------------------------------
  # 1. UMAP by Batch
  # ---------------------------------------------------------------------------
  if (batch_var %in% available_vars) {
    log_message("  Plotting UMAP by batch...")
    
    n_batches <- length(unique(obj[[batch_var, drop = TRUE]]))
    
    p_batch <- DimPlot(obj, 
                       reduction = "umap", 
                       group.by = batch_var, 
                       pt.size = pt_size,
                       raster = use_raster,
                       raster.dpi = raster_dpi) +
      labs(title = paste0("Batch: ", batch_var),
           subtitle = paste0(format(ncol(obj), big.mark = ","), " cells | ",
                             n_batches, " batches")) +
      scale_color_viridis(option = "H", discrete = TRUE) +
      guides(color = "none") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_batch.png")),
           p_batch, width = 10, height = 8, dpi = 300, bg = "white")
  }
  
  # ---------------------------------------------------------------------------
  # 2. UMAP by Cell Type L1
  # ---------------------------------------------------------------------------
  if (labels_l1 %in% available_vars) {
    log_message("  Plotting UMAP by cell type L1...")
    
    n_types <- length(unique(obj[[labels_l1, drop = TRUE]]))
    colors_l1 <- colorRampPalette(brewer.pal(12, "Paired"))(n_types)
    
    p_l1 <- DimPlot(obj, 
                    reduction = "umap", 
                    group.by = labels_l1, 
                    pt.size = pt_size,
                    label = TRUE,
                    label.size = 3,
                    repel = TRUE,
                    raster = use_raster,
                    raster.dpi = raster_dpi) +
      labs(title = "Cell Type (L1)") +
      scale_color_manual(values = colors_l1) +
      guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_celltype_l1.png")),
           p_l1, width = 12, height = 8, dpi = 300, bg = "white")
  }
  
  # ---------------------------------------------------------------------------
  # 3. UMAP by Cell Type L2
  # ---------------------------------------------------------------------------
  if (labels_l2 %in% available_vars) {
    log_message("  Plotting UMAP by cell type L2...")
    
    n_types <- length(unique(obj[[labels_l2, drop = TRUE]]))
    colors_l2 <- colorRampPalette(brewer.pal(12, "Set3"))(n_types)
    
    p_l2 <- DimPlot(obj, 
                    reduction = "umap", 
                    group.by = labels_l2, 
                    pt.size = pt_size,
                    label = TRUE,
                    label.size = 2,
                    repel = TRUE,
                    raster = use_raster,
                    raster.dpi = raster_dpi) +
      labs(title = "Cell Type (L2)") +
      scale_color_manual(values = colors_l2) +
      guides(color = guide_legend(override.aes = list(size = 2), ncol = 2)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            legend.text = element_text(size = 6))
    
    ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_celltype_l2.png")),
           p_l2, width = 14, height = 8, dpi = 300, bg = "white")
  }
  
  # ---------------------------------------------------------------------------
  # 4. UMAP by Cluster
  # ---------------------------------------------------------------------------
  if ("leiden.cluster" %in% available_vars) {
    log_message("  Plotting UMAP by cluster...")
    
    n_clusters <- length(unique(obj$leiden.cluster))
    
    p_cluster <- DimPlot(obj, 
                         reduction = "umap", 
                         group.by = "leiden.cluster", 
                         pt.size = pt_size,
                         label = TRUE,
                         label.size = 4,
                         raster = use_raster,
                         raster.dpi = raster_dpi) +
      labs(title = "Leiden Clusters",
           subtitle = paste0(n_clusters, " clusters")) +
      scale_color_viridis_d(option = "turbo") +
      guides(color = "none") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_clusters.png")),
           p_cluster, width = 10, height = 8, dpi = 300, bg = "white")
  }
  
  # ---------------------------------------------------------------------------
  # 5. UMAP by Tissue
  # ---------------------------------------------------------------------------
  if ("tissue" %in% available_vars) {
    log_message("  Plotting UMAP by tissue...")
    
    tissue_colors <- c("Tumor" = "#E41A1C", "Normal" = "#377EB8", 
                       "Blood" = "#4DAF4A", "LN" = "#984EA3", 
                       "Met" = "#FF7F00", "Juxta" = "#FFFF33")
    
    p_tissue <- DimPlot(obj, 
                        reduction = "umap", 
                        group.by = "tissue", 
                        pt.size = pt_size,
                        raster = use_raster,
                        raster.dpi = raster_dpi) +
      labs(title = "Tissue Type") +
      scale_color_manual(values = tissue_colors, na.value = "grey80") +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_tissue.png")),
           p_tissue, width = 10, height = 8, dpi = 300, bg = "white")
  }
  
  # ---------------------------------------------------------------------------
  # 6. UMAP by Clonal Expansion
  # ---------------------------------------------------------------------------
  if ("cloneSize" %in% available_vars) {
    log_message("  Plotting UMAP by clonal expansion...")
    
    p_clone <- DimPlot(obj, 
                       reduction = "umap", 
                       group.by = "cloneSize", 
                       pt.size = pt_size,
                       raster = use_raster,
                       raster.dpi = raster_dpi) +
      labs(title = "Clonal Expansion") +
      scale_color_viridis(option = "B", discrete = TRUE, na.value = "grey") +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_clonalExpansion.png")),
           p_clone, width = 11, height = 8, dpi = 300, bg = "white")
  }
  
  # ---------------------------------------------------------------------------
  # 7. Feature Plots (T cell markers)
  # ---------------------------------------------------------------------------
  log_message("  Plotting T cell markers...")
    
  markers <- c("CD3E", "CD4", "CD8A", "CD8B", 
               "FOXP3", "CCR7", "TCF7", "SELL",
               "GZMB", "PRF1", "IFNG", "PDCD1",
               "HAVCR2", "LAG3", "TIGIT", "TOX")
    
  available_markers <- intersect(markers, rownames(obj))
    
  if (length(available_markers) >= 4) {
      for (i in seq(1, min(length(available_markers), 12), by = 4)) {
        batch_markers <- available_markers[i:min(i+3, length(available_markers))]
        
        p_feat <- FeaturePlot(obj, 
                              features = batch_markers, 
                              reduction = "umap",
                              pt.size = pt_size * 0.5,
                              ncol = 2,
                              raster = use_raster,
                              raster.dpi = raster_dpi) &
          scale_color_viridis_c() &
          theme_minimal()
        
        batch_num <- ceiling(i / 4)
        ggsave(file.path(fig_dir, paste0(output_prefix, "_umap_markers_", batch_num, ".png")),
               p_feat, width = 10, height = 10, dpi = 300, bg = "white")
      }
    }
  
  log_message("--- Plots saved to: ", fig_dir, " ---")
}