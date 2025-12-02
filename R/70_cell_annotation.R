# 50_cell_annotation.R
# ==============================================================================
# Reproducible T Cell Annotation System
# ==============================================================================
# 
# This module provides a hierarchical, marker-based annotation system for T cells
# similar to scGate, combined with consensus voting from reference annotations.
#
# Key Features:
# - Hierarchical gating based on canonical T cell markers
# - Integration with Azimuth (predicted.celltype), HPCA, and Monaco labels
# - Cluster-level consensus annotation
# - Comprehensive visualization functions
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(pheatmap)
})

# ==============================================================================
# 1. MARKER LOADING AND PARSING
# ==============================================================================

#' Load T cell marker definitions from YAML
#'
#' @param yaml_path Path to the tcell_markers.yaml file
#' @return List containing marker definitions and parameters
load_tcell_markers <- function(yaml_path = "tcell_markers.yaml") {
  if (!file.exists(yaml_path)) {
    stop("Marker definition file not found: ", yaml_path)
  }
  yaml::read_yaml(yaml_path)
}

#' Get all markers from a category
#'
#' @param markers_config Loaded marker configuration
#' @param category Category name (e.g., "lineages", "cd4_subtypes")
#' @return Named list of marker sets per cell type
get_markers_by_category <- function(markers_config, category) {
  if (!category %in% names(markers_config)) {
    stop("Category not found: ", category)
  }
  markers_config[[category]]
}

#' Flatten all unique markers from config
#'
#' @param markers_config Loaded marker configuration
#' @return Character vector of all unique markers
get_all_markers <- function(markers_config) {
  all_markers <- c()
  
  categories <- c("lineages", "cd4_subtypes", "cd8_subtypes", 
                  "other_states", "special", "nk_subtypes")
  
  for (cat in categories) {
    if (cat %in% names(markers_config)) {
      for (cell_type in names(markers_config[[cat]])) {
        ct_info <- markers_config[[cat]][[cell_type]]
        all_markers <- c(all_markers, 
                        ct_info$positive, 
                        ct_info$negative, 
                        ct_info$optional)
      }
    }
  }
  
  unique(all_markers[!is.na(all_markers)])
}

# ==============================================================================
# 2. MARKER EXPRESSION SCORING
# ==============================================================================

#' Calculate marker module scores for each cell
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param assay Assay to use for expression
#' @return Seurat object with module scores added
calculate_marker_scores <- function(obj, markers_config, assay = "RNA") {
  
  log_message("Calculating marker module scores...")
  
  DefaultAssay(obj) <- assay
  
  # Get all markers and check availability

  all_markers <- get_all_markers(markers_config)
  available_markers <- intersect(all_markers, rownames(obj))
  missing_markers <- setdiff(all_markers, available_markers)
  
  if (length(missing_markers) > 0) {
    log_message("Note: ", length(missing_markers), " markers not found in data")
    log_message("Missing: ", paste(head(missing_markers, 10), collapse = ", "), 
                if(length(missing_markers) > 10) "..." else "")
  }
  
  # Calculate scores for each cell type
  categories <- c("lineages", "cd4_subtypes", "cd8_subtypes", 
                  "other_states", "nk_subtypes")
  
  for (cat in categories) {
    if (!cat %in% names(markers_config)) next
    
    for (cell_type in names(markers_config[[cat]])) {
      ct_info <- markers_config[[cat]][[cell_type]]
      
      # Get positive markers that exist in data
      pos_markers <- intersect(ct_info$positive, available_markers)
      
      if (length(pos_markers) >= 1) {
        score_name <- paste0("Score_", cell_type)
        
        # Use AddModuleScore for positive markers
        tryCatch({
          obj <- AddModuleScore(
            obj,
            features = list(pos_markers),
            name = score_name,
            assay = assay,
            search = FALSE
          )
          # Rename the score column (AddModuleScore adds "1" suffix)
          old_name <- paste0(score_name, "1")
          if (old_name %in% colnames(obj@meta.data)) {
            obj@meta.data[[score_name]] <- obj@meta.data[[old_name]]
            obj@meta.data[[old_name]] <- NULL
          }
        }, error = function(e) {
          log_message("Warning: Could not calculate score for ", cell_type)
        })
      }
    }
  }
  
  obj
}

#' Calculate expression fraction per cluster for markers
#'
#' @param obj Seurat object
#' @param markers Character vector of marker genes
#' @param cluster_col Column name for clusters
#' @param threshold Expression threshold (default 0)
#' @return Matrix of expression fractions (clusters x markers)
calculate_cluster_marker_expression <- function(obj, markers, cluster_col = "leiden.cluster", 
                                                  threshold = 0) {
  
  DefaultAssay(obj) <- "RNA"
  
  # Get available markers
  available_markers <- intersect(markers, rownames(obj))
  
  if (length(available_markers) == 0) {
    warning("No markers found in data")
    return(NULL)
  }
  
  # Get expression matrix
  expr_mat <- GetAssayData(obj, layer = "data")[available_markers, , drop = FALSE]
  clusters <- obj@meta.data[[cluster_col]]
  
  # Calculate fraction expressing per cluster
  cluster_levels <- sort(unique(clusters))
  result <- matrix(0, nrow = length(cluster_levels), ncol = length(available_markers),
                   dimnames = list(cluster_levels, available_markers))
  
  for (cl in cluster_levels) {
    cells <- which(clusters == cl)
    if (length(cells) > 0) {
      for (marker in available_markers) {
        expr <- expr_mat[marker, cells]
        result[as.character(cl), marker] <- mean(expr > threshold)
      }
    }
  }
  
  result
}

#' Calculate mean expression per cluster for markers
#'
#' @param obj Seurat object
#' @param markers Character vector of marker genes
#' @param cluster_col Column name for clusters
#' @return Matrix of mean expressions (clusters x markers)
calculate_cluster_mean_expression <- function(obj, markers, cluster_col = "leiden.cluster") {
  
  DefaultAssay(obj) <- "RNA"
  
  available_markers <- intersect(markers, rownames(obj))
  
  if (length(available_markers) == 0) {
    warning("No markers found in data")
    return(NULL)
  }
  
  expr_mat <- GetAssayData(obj, layer = "data")[available_markers, , drop = FALSE]
  clusters <- obj@meta.data[[cluster_col]]
  
  cluster_levels <- sort(unique(clusters))
  result <- matrix(0, nrow = length(cluster_levels), ncol = length(available_markers),
                   dimnames = list(cluster_levels, available_markers))
  
  for (cl in cluster_levels) {
    cells <- which(clusters == cl)
    if (length(cells) > 0) {
      result[as.character(cl), ] <- rowMeans(expr_mat[, cells, drop = FALSE])
    }
  }
  
  result
}

# ==============================================================================
# 3. HIERARCHICAL GATING
# ==============================================================================

#' Apply scGate-style hierarchical gating to cells
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param assay Assay to use
#' @return Seurat object with gating annotations
apply_hierarchical_gating <- function(obj, markers_config, assay = "RNA") {
  
  log_message("Applying hierarchical gating...")
  
  DefaultAssay(obj) <- assay
  scoring_params <- markers_config$scoring
  
  # Initialize annotation columns
  obj$gate_lineage <- NA_character_
  obj$gate_subtype <- NA_character_
  obj$gate_state <- NA_character_
  
  # Get expression matrix
  all_markers <- get_all_markers(markers_config)
  available_markers <- intersect(all_markers, rownames(obj))
  expr_mat <- GetAssayData(obj, layer = "data")[available_markers, , drop = FALSE]
  
  # Calculate expression percentiles for thresholding
  expr_percentiles <- apply(expr_mat, 1, function(x) quantile(x[x > 0], 0.5, na.rm = TRUE))
  expr_percentiles[is.na(expr_percentiles)] <- 0
  
  # Helper function to score a cell type
  score_cell_type <- function(cell_idx, ct_info) {
    pos_markers <- intersect(ct_info$positive, available_markers)
    neg_markers <- intersect(ct_info$negative, available_markers)
    
    if (length(pos_markers) == 0) return(0)
    
    # Positive score: fraction of positive markers expressed
    pos_expr <- sapply(pos_markers, function(m) {
      expr_mat[m, cell_idx] > expr_percentiles[m]
    })
    pos_score <- mean(pos_expr)
    
    # Negative score: penalty for expressing negative markers
    neg_score <- 0
    if (length(neg_markers) > 0) {
      neg_expr <- sapply(neg_markers, function(m) {
        expr_mat[m, cell_idx] > expr_percentiles[m]
      })
      neg_score <- mean(neg_expr)
    }
    
    # Combined score
    final_score <- pos_score * scoring_params$positive_weight - 
                   neg_score * scoring_params$negative_weight
    
    # Check minimum positive fraction requirement
    if (pos_score < scoring_params$min_positive_fraction) {
      final_score <- final_score * 0.5  # Penalize
    }
    
    final_score
  }
  
  # --- Level 1: Lineage Assignment ---
  lineages <- markers_config$lineages
  n_cells <- ncol(obj)
  
  log_message("  Level 1: Assigning lineages (", n_cells, " cells)...")
  
  lineage_scores <- sapply(names(lineages), function(lin) {
    sapply(1:n_cells, function(i) score_cell_type(i, lineages[[lin]]))
  })
  
  # Assign lineage based on highest score
  obj$gate_lineage <- apply(lineage_scores, 1, function(row) {
    if (max(row) > 0) names(row)[which.max(row)] else "Unknown"
  })
  obj$gate_lineage_score <- apply(lineage_scores, 1, max)
  
  # --- Level 2: Subtype Assignment ---
  log_message("  Level 2: Assigning subtypes...")
  
  # CD4 subtypes
  cd4_cells <- which(obj$gate_lineage == "CD4_T")
  if (length(cd4_cells) > 0) {
    cd4_subtypes <- markers_config$cd4_subtypes
    subtype_scores <- sapply(names(cd4_subtypes), function(st) {
      sapply(cd4_cells, function(i) score_cell_type(i, cd4_subtypes[[st]]))
    })
    
    if (!is.matrix(subtype_scores)) {
      subtype_scores <- matrix(subtype_scores, ncol = length(cd4_subtypes),
                               dimnames = list(NULL, names(cd4_subtypes)))
    }
    
    obj$gate_subtype[cd4_cells] <- apply(subtype_scores, 1, function(row) {
      if (max(row) > 0) names(row)[which.max(row)] else "CD4_unassigned"
    })
  }
  
  # CD8 subtypes
  cd8_cells <- which(obj$gate_lineage == "CD8_T")
  if (length(cd8_cells) > 0) {
    cd8_subtypes <- markers_config$cd8_subtypes
    subtype_scores <- sapply(names(cd8_subtypes), function(st) {
      sapply(cd8_cells, function(i) score_cell_type(i, cd8_subtypes[[st]]))
    })
    
    if (!is.matrix(subtype_scores)) {
      subtype_scores <- matrix(subtype_scores, ncol = length(cd8_subtypes),
                               dimnames = list(NULL, names(cd8_subtypes)))
    }
    
    obj$gate_subtype[cd8_cells] <- apply(subtype_scores, 1, function(row) {
      if (max(row) > 0) names(row)[which.max(row)] else "CD8_unassigned"
    })
  }
  
  # NK subtypes
  nk_cells <- which(obj$gate_lineage == "NK")
  if (length(nk_cells) > 0 && "nk_subtypes" %in% names(markers_config)) {
    nk_subtypes <- markers_config$nk_subtypes
    subtype_scores <- sapply(names(nk_subtypes), function(st) {
      sapply(nk_cells, function(i) score_cell_type(i, nk_subtypes[[st]]))
    })
    
    if (!is.matrix(subtype_scores)) {
      subtype_scores <- matrix(subtype_scores, ncol = length(nk_subtypes),
                               dimnames = list(NULL, names(nk_subtypes)))
    }
    
    obj$gate_subtype[nk_cells] <- apply(subtype_scores, 1, function(row) {
      if (max(row) > 0) names(row)[which.max(row)] else "NK_unassigned"
    })
  }
  
  # Other lineages get their lineage as subtype
  other_cells <- which(obj$gate_lineage %in% c("gdT", "MAIT", "NKT", "dnT"))
  obj$gate_subtype[other_cells] <- obj$gate_lineage[other_cells]
  
  # --- Level 3: State Assignment ---
  log_message("  Level 3: Assigning cellular states...")
  
  other_states <- markers_config$other_states
  special_states <- markers_config$special
  all_states <- c(other_states, special_states)
  
  state_scores <- sapply(names(all_states), function(st) {
    sapply(1:n_cells, function(i) score_cell_type(i, all_states[[st]]))
  })
  
  # Assign state if above threshold
  obj$gate_state <- apply(state_scores, 1, function(row) {
    if (max(row) > 0.3) names(row)[which.max(row)] else "None"
  })
  
  log_message("  Gating complete.")
  obj
}

# ==============================================================================
# 4. CLUSTER-LEVEL CONSENSUS ANNOTATION
# ==============================================================================

#' Get majority vote annotation per cluster
#'
#' @param obj Seurat object
#' @param annotation_col Column with cell annotations
#' @param cluster_col Column with cluster IDs
#' @param min_fraction Minimum fraction for majority (default 0.3)
#' @return Data frame with cluster annotations
get_cluster_majority <- function(obj, annotation_col, cluster_col = "leiden.cluster",
                                  min_fraction = 0.3) {
  
  df <- obj@meta.data %>%
    select(all_of(c(cluster_col, annotation_col))) %>%
    filter(!is.na(.data[[annotation_col]])) %>%
    group_by(.data[[cluster_col]], .data[[annotation_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(
      total = sum(count),
      fraction = count / total
    ) %>%
    slice_max(count, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Mark low-confidence assignments
  df$confident <- df$fraction >= min_fraction
  
  df
}

#' Generate consensus cluster annotations from multiple sources
#'
#' @param obj Seurat object
#' @param cluster_col Column with cluster IDs
#' @param reference_cols Reference annotation columns to use
#' @param markers_config Loaded marker configuration (for gating results)
#' @return Data frame with consensus annotations per cluster
generate_cluster_consensus <- function(obj, 
                                        cluster_col = "leiden.cluster",
                                        reference_cols = c("predicted.celltype.l2", 
                                                           "Monaco.labels",
                                                           "gate_subtype"),
                                        markers_config = NULL) {
  
  log_message("Generating cluster consensus annotations...")
  
  clusters <- sort(unique(obj@meta.data[[cluster_col]]))
  
  # Get majority vote from each reference
  majorities <- list()
  for (ref_col in reference_cols) {
    if (ref_col %in% colnames(obj@meta.data)) {
      majorities[[ref_col]] <- get_cluster_majority(obj, ref_col, cluster_col)
    }
  }
  
  # Build consensus table
  consensus_df <- data.frame(cluster = clusters)
  
  for (ref_col in names(majorities)) {
    maj <- majorities[[ref_col]]
    consensus_df[[paste0(ref_col, "_label")]] <- 
      maj[[ref_col]][match(consensus_df$cluster, maj[[cluster_col]])]
    consensus_df[[paste0(ref_col, "_fraction")]] <- 
      maj$fraction[match(consensus_df$cluster, maj[[cluster_col]])]
  }
  
  # Calculate cluster sizes
  cluster_sizes <- table(obj@meta.data[[cluster_col]])
  consensus_df$n_cells <- as.numeric(cluster_sizes[as.character(consensus_df$cluster)])
  
  consensus_df
}

#' Assign final cluster annotation based on decision rules
#'
#' @param consensus_df Output from generate_cluster_consensus
#' @param priority_order Priority of reference annotations
#' @return Data frame with final annotations
assign_final_cluster_annotation <- function(consensus_df,
                                             priority_order = c("gate_subtype",
                                                                "predicted.celltype.l2",
                                                                "Monaco.labels")) {
  
  # Create mapping for consistent naming
  type_mapping <- list(
    # CD4 subtypes
    "Th1" = "Th1",
    "Th1 cells" = "Th1",
    "Th2" = "Th2", 
    "Th2 cells" = "Th2",
    "Th17" = "Th17",
    "Th17 cells" = "Th17",
    "Th1/Th17 cells" = "Th1/Th17",
    "Treg" = "Treg",
    "Treg_naive" = "Treg_naive",
    "Treg_memory" = "Treg_memory",
    "T regulatory cells" = "Treg",
    "Tfh" = "Tfh",
    "Follicular helper T cells" = "Tfh",
    
    # CD4 memory/naive
    "CD4 Naive" = "CD4_naive",
    "Naive CD4 T cells" = "CD4_naive",
    "CD4 TCM" = "CD4_Tcm",
    "CD4 TEM" = "CD4_Tem",
    "CD4 CTL" = "CD4_CTL",
    
    # CD8 subtypes
    "CD8 Naive" = "CD8_naive",
    "Naive CD8 T cells" = "CD8_naive",
    "CD8 TCM" = "CD8_Tcm",
    "Central memory CD8 T cells" = "CD8_Tcm",
    "CD8 TEM" = "CD8_Tem",
    "Effector memory CD8 T cells" = "CD8_Tem",
    "CD8_Temra" = "CD8_Temra",
    "Terminal effector CD8 T cells" = "CD8_Temra",
    "CD8_exhausted" = "CD8_exhausted",
    "CD8_CTL" = "CD8_CTL",
    
    # Other
    "MAIT" = "MAIT",
    "MAIT cells" = "MAIT",
    "gdT" = "gdT",
    "Vd2 gd T cells" = "gdT_Vd2",
    "Non-Vd2 gd T cells" = "gdT_nonVd2",
    "dnT" = "dnT",
    "NK" = "NK",
    "NK_CD56bright" = "NK_CD56bright",
    "NK_CD56dim" = "NK_CD56dim",
    "Natural killer cells" = "NK"
  )
  
  # Determine final annotation per cluster
  final_annotations <- sapply(1:nrow(consensus_df), function(i) {
    row <- consensus_df[i, ]
    
    best_label <- NA
    best_confidence <- 0
    
    for (ref in priority_order) {
      label_col <- paste0(ref, "_label")
      frac_col <- paste0(ref, "_fraction")
      
      if (label_col %in% names(row) && !is.na(row[[label_col]])) {
        frac <- ifelse(!is.na(row[[frac_col]]), row[[frac_col]], 0)
        
        # Use gate_subtype with priority if high confidence
        if (ref == "gate_subtype" && frac > 0.5) {
          best_label <- row[[label_col]]
          best_confidence <- frac
          break
        }
        
        # Otherwise take highest confidence
        if (frac > best_confidence) {
          best_label <- row[[label_col]]
          best_confidence <- frac
        }
      }
    }
    
    # Map to canonical name
    if (!is.na(best_label) && best_label %in% names(type_mapping)) {
      best_label <- type_mapping[[best_label]]
    }
    
    best_label
  })
  
  consensus_df$final_annotation <- final_annotations
  consensus_df$annotation_confidence <- sapply(1:nrow(consensus_df), function(i) {
    fracs <- consensus_df[i, grep("_fraction$", names(consensus_df))]
    mean(as.numeric(fracs), na.rm = TRUE)
  })
  
  consensus_df
}

#' Apply cluster annotations to Seurat object
#'
#' @param obj Seurat object
#' @param consensus_df Consensus data frame with final annotations
#' @param cluster_col Column with cluster IDs
#' @return Seurat object with new annotation column
apply_cluster_annotations <- function(obj, consensus_df, cluster_col = "leiden.cluster") {
  
  # Create mapping
  annotation_map <- setNames(consensus_df$final_annotation, 
                              as.character(consensus_df$cluster))
  
  # Apply to object

  obj$cluster_annotation <- annotation_map[as.character(obj@meta.data[[cluster_col]])]
  
  # Also add confidence
  confidence_map <- setNames(consensus_df$annotation_confidence,
                              as.character(consensus_df$cluster))
  obj$annotation_confidence <- confidence_map[as.character(obj@meta.data[[cluster_col]])]
  
  obj
}

# ==============================================================================
# 5. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Plot canonical T cell markers as DotPlot
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param group.by Grouping variable
#' @return ggplot object
plot_canonical_markers <- function(obj, markers_config, group.by = "cluster_annotation") {
  
  # Define key markers for each major type
  key_markers <- list(
    "Pan-T" = c("CD3D", "CD3E"),
    "CD4" = c("CD4"),
    "CD8" = c("CD8A", "CD8B"),
    "Naive" = c("CCR7", "SELL", "LEF1", "TCF7"),
    "Memory" = c("IL7R", "CD27"),
    "Effector" = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7"),
    "Th1" = c("TBX21", "IFNG", "CXCR3"),
    "Th2" = c("GATA3", "IL4", "CCR4"),
    "Th17" = c("RORC", "IL17A", "CCR6"),
    "Tfh" = c("BCL6", "CXCR5", "PDCD1", "ICOS"),
    "Treg" = c("FOXP3", "IL2RA", "CTLA4"),
    "Exhausted" = c("PDCD1", "LAG3", "HAVCR2", "TOX"),
    "NK" = c("NCAM1", "KLRD1", "KLRF1"),
    "MAIT" = c("KLRB1", "SLC4A10"),
    "gdT" = c("TRDC", "TRGC1"),
    "Proliferation" = c("MKI67", "TOP2A")
  )
  
  # Get available markers
  all_markers <- unlist(key_markers)
  available <- intersect(all_markers, rownames(obj))
  
  # Create ordered marker list
  marker_order <- c()
  for (cat in names(key_markers)) {
    avail_in_cat <- intersect(key_markers[[cat]], available)
    marker_order <- c(marker_order, avail_in_cat)
  }
  
  # Create DotPlot
  p <- DotPlot(obj, features = marker_order, group.by = group.by) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_viridis_c() +
    labs(title = "Canonical T Cell Markers")
  
  p
}

#' Plot cluster annotation comparison heatmap
#'
#' @param consensus_df Consensus data frame
#' @return Heatmap object
plot_annotation_comparison <- function(consensus_df) {
  
  # Extract label columns
  label_cols <- grep("_label$", names(consensus_df), value = TRUE)
  
  if (length(label_cols) < 2) {
    warning("Need at least 2 annotation sources for comparison")
    return(NULL)
  }
  
  # Create factor matrix for heatmap
  labels_df <- consensus_df[, label_cols]
  rownames(labels_df) <- paste0("Cluster_", consensus_df$cluster)
  
  # Add final annotation
  labels_df$final <- consensus_df$final_annotation
  
  # Convert to numeric for heatmap (all unique labels)
  all_labels <- unique(unlist(labels_df))
  label_numeric <- lapply(labels_df, function(x) match(x, all_labels))
  label_mat <- do.call(cbind, label_numeric)
  rownames(label_mat) <- rownames(labels_df)
  
  # Create color annotation
  n_labels <- length(all_labels)
  label_colors <- setNames(
    colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                       "#FF7F00", "#FFFF33", "#A65628", "#F781BF"))(n_labels),
    all_labels
  )
  
  # Plot
  pheatmap(
    label_mat,
    color = label_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Annotation Source Comparison",
    legend = FALSE,
    annotation_names_row = TRUE
  )
}

#' Create UMAP with annotations
#'
#' @param obj Seurat object
#' @param group.by Grouping variable(s)
#' @param reduction Reduction to use
#' @return ggplot object
plot_annotation_umap <- function(obj, group.by = "cluster_annotation", 
                                  reduction = "umap") {
  
  plots <- list()
  
  for (gb in group.by) {
    if (gb %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = reduction, group.by = gb, 
                   label = TRUE, label.size = 3, repel = TRUE) +
        theme(legend.position = "right") +
        labs(title = gb)
      plots[[gb]] <- p
    }
  }
  
  if (length(plots) > 0) {
    wrap_plots(plots, ncol = 2)
  } else {
    warning("No valid grouping variables found")
    NULL
  }
}

#' Plot marker expression heatmap per cluster
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param cluster_col Cluster column
#' @param scale Scale expression values
#' @return ComplexHeatmap object
plot_marker_heatmap <- function(obj, markers_config, 
                                 cluster_col = "cluster_annotation",
                                 scale = TRUE) {
  
  # Get key markers
  key_markers <- c(
    # Lineage
    "CD3D", "CD4", "CD8A", "CD8B", "NCAM1",
    # Naive/Memory
    "CCR7", "SELL", "LEF1", "TCF7", "IL7R",
    # Effector
    "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7",
    # CD4 subtypes
    "TBX21", "IFNG", "GATA3", "RORC", "IL17A", "BCL6", "CXCR5", "FOXP3", "IL2RA",
    # Exhaustion
    "PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT",
    # Other
    "MKI67", "KLRB1", "TRDC"
  )
  
  available <- intersect(key_markers, rownames(obj))
  
  # Calculate mean expression per cluster
  expr_mat <- calculate_cluster_mean_expression(obj, available, cluster_col)
  
  if (is.null(expr_mat)) return(NULL)
  
  # Scale if requested
  if (scale) {
    expr_mat <- t(scale(t(expr_mat)))
    expr_mat[expr_mat > 2] <- 2
    expr_mat[expr_mat < -2] <- -2
  }
  
  # Create heatmap
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  Heatmap(
    t(expr_mat),
    name = "Scaled\nExpression",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "Clusters",
    row_title = "Markers"
  )
}

# ==============================================================================
# 6. MAIN ANNOTATION PIPELINE
# ==============================================================================

#' Run complete T cell annotation pipeline
#'
#' @param obj Seurat object
#' @param markers_yaml Path to marker YAML file
#' @param cluster_col Cluster column name
#' @param output_dir Output directory for results
#' @return List with annotated object and summary
run_tcell_annotation <- function(obj, 
                                  markers_yaml = "tcell_markers.yaml",
                                  cluster_col = "leiden.cluster",
                                  output_dir = NULL) {
  
  log_message("=== Starting T Cell Annotation Pipeline ===")
  
  # 1. Load marker definitions
  markers_config <- load_tcell_markers(markers_yaml)
  log_message("Loaded marker definitions")
  
  # 2. Calculate marker module scores
  obj <- calculate_marker_scores(obj, markers_config)
  
  # 3. Apply hierarchical gating
  obj <- apply_hierarchical_gating(obj, markers_config)
  
  # 4. Generate cluster consensus
  consensus_df <- generate_cluster_consensus(
    obj, 
    cluster_col = cluster_col,
    reference_cols = c("predicted.celltype.l2", "Monaco.labels", "gate_subtype"),
    markers_config = markers_config
  )
  
  # 5. Assign final annotations
  consensus_df <- assign_final_cluster_annotation(consensus_df)
  
  # 6. Apply to object
  obj <- apply_cluster_annotations(obj, consensus_df, cluster_col)
  
  # 7. Generate summary
  annotation_summary <- table(obj$cluster_annotation) %>%
    as.data.frame() %>%
    arrange(desc(Freq)) %>%
    setNames(c("Annotation", "Cell_Count"))
  
  log_message("Annotation complete!")
  log_message("Summary of annotations:")
  print(annotation_summary)
  
  # 8. Save outputs if directory provided
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save consensus table
    write.csv(consensus_df, file.path(output_dir, "cluster_annotations.csv"), 
              row.names = FALSE)
    
    # Save summary
    write.csv(annotation_summary, file.path(output_dir, "annotation_summary.csv"),
              row.names = FALSE)
    
    log_message("Results saved to: ", output_dir)
  }
  
  list(
    object = obj,
    consensus = consensus_df,
    summary = annotation_summary,
    markers_config = markers_config
  )
}

# ==============================================================================
# 7. QUALITY CONTROL AND VALIDATION
# ==============================================================================

#' Validate annotations against known markers
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param annotation_col Annotation column to validate
#' @return Data frame with validation scores
validate_annotations <- function(obj, markers_config, 
                                  annotation_col = "cluster_annotation") {
  
  annotations <- unique(obj@meta.data[[annotation_col]])
  annotations <- annotations[!is.na(annotations)]
  
  validation_results <- data.frame()
  
  for (annot in annotations) {
    cells <- which(obj@meta.data[[annotation_col]] == annot)
    
    if (length(cells) == 0) next
    
    # Find expected markers for this annotation
    expected_pos <- c()
    expected_neg <- c()
    
    # Check all categories
    for (cat in c("lineages", "cd4_subtypes", "cd8_subtypes", "nk_subtypes")) {
      if (cat %in% names(markers_config)) {
        for (ct in names(markers_config[[cat]])) {
          if (grepl(ct, annot, ignore.case = TRUE) || 
              grepl(annot, ct, ignore.case = TRUE)) {
            expected_pos <- c(expected_pos, markers_config[[cat]][[ct]]$positive)
            expected_neg <- c(expected_neg, markers_config[[cat]][[ct]]$negative)
          }
        }
      }
    }
    
    expected_pos <- intersect(unique(expected_pos), rownames(obj))
    expected_neg <- intersect(unique(expected_neg), rownames(obj))
    
    # Calculate expression in these cells
    if (length(expected_pos) > 0) {
      expr <- GetAssayData(obj, layer = "data")[expected_pos, cells, drop = FALSE]
      pos_frac <- mean(rowMeans(expr > 0))
    } else {
      pos_frac <- NA
    }
    
    if (length(expected_neg) > 0) {
      expr <- GetAssayData(obj, layer = "data")[expected_neg, cells, drop = FALSE]
      neg_frac <- mean(rowMeans(expr > 0))
    } else {
      neg_frac <- NA
    }
    
    validation_results <- rbind(validation_results, data.frame(
      annotation = annot,
      n_cells = length(cells),
      n_pos_markers = length(expected_pos),
      n_neg_markers = length(expected_neg),
      pos_marker_expr_frac = pos_frac,
      neg_marker_expr_frac = neg_frac,
      validation_score = ifelse(!is.na(pos_frac) & !is.na(neg_frac),
                                pos_frac - neg_frac,
                                ifelse(!is.na(pos_frac), pos_frac, 0))
    ))
  }
  
  validation_results %>% arrange(desc(validation_score))
}

#' Generate annotation quality report
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param output_file Output file path
#' @return Invisible NULL
generate_annotation_report <- function(obj, markers_config, 
                                        output_file = "annotation_report.html") {
  
  log_message("Generating annotation quality report...")
  
  # Validation results
  validation <- validate_annotations(obj, markers_config)
  
  # Create simple text report for now
  report_text <- c(
    "# T Cell Annotation Quality Report",
    "",
    paste("Date:", Sys.Date()),
    paste("Total cells:", ncol(obj)),
    paste("Total clusters:", length(unique(obj$leiden.cluster))),
    "",
    "## Annotation Summary",
    "",
    capture.output(print(table(obj$cluster_annotation))),
    "",
    "## Validation Scores",
    "",
    capture.output(print(validation))
  )
  
  writeLines(report_text, gsub("\\.html$", ".txt", output_file))
  log_message("Report saved to: ", gsub("\\.html$", ".txt", output_file))
  
  invisible(NULL)
}

# ==============================================================================
# 8. UTILITY FUNCTIONS
# ==============================================================================

#' Log message with timestamp
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
}

#' Export annotations for external tools
#'
#' @param obj Seurat object
#' @param output_file Output CSV file
#' @param cols Columns to export
export_annotations <- function(obj, output_file, 
                                cols = c("leiden.cluster", "cluster_annotation",
                                         "gate_lineage", "gate_subtype",
                                         "predicted.celltype.l2", "Monaco.labels")) {
  
  available_cols <- intersect(cols, colnames(obj@meta.data))
  
  export_df <- obj@meta.data[, available_cols, drop = FALSE]
  export_df$cell_barcode <- rownames(obj@meta.data)
  
  write.csv(export_df, output_file, row.names = FALSE)
  log_message("Annotations exported to: ", output_file)
}
