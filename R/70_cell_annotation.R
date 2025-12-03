# 70_cell_annotation.R
# ==============================================================================
# Reproducible T Cell Annotation System 
# ==============================================================================


suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(RColorBrewer)
  library(scales)
  library(grid)
  library(gridExtra)
})

# Try to load ComplexHeatmap, fall back to pheatmap if not available
use_complex_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
if (use_complex_heatmap) {
  library(ComplexHeatmap)
  library(circlize)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Log message with timestamp
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
}

#' Get publication-quality color palette for T cell annotations
get_tcell_palette <- function() {
  c(
    # CD4 subtypes - Red/Orange spectrum
    "CD4_T" = "#E41A1C",
    "CD4_naive" = "#FFC0CB",
    "CD4_Tcm" = "#FF6B6B",
    "CD4_Tem" = "#DC143C",
    "CD4_Temra" = "#8B0000",
    "CD4_CTL" = "#B22222",
    "CD4_proliferating" = "#FF4500",
    
    # Th subtypes - Warm colors
    "Th1" = "#FF8C00",
    "Th2" = "#FFA500",
    "Th17" = "#FF7F00",
    "Th1_Th17" = "#FF6347",
    "Tfh" = "#9370DB",
    
    # Treg - Gold/Yellow
    "Treg" = "#FFD700",
    "Treg_naive" = "#FFEC8B",
    "Treg_memory" = "#DAA520",
    "Treg_activated" = "#B8860B",
    
    # CD8 subtypes - Blue spectrum
    "CD8_T" = "#4169E1",
    "CD8_naive" = "#87CEEB",
    "CD8_Tcm" = "#6495ED",
    "CD8_Tem" = "#4682B4",
    "CD8_Temra" = "#000080",
    "CD8_exhausted" = "#191970",
    "CD8_CTL" = "#0000CD",
    "CD8_proliferating" = "#1E90FF",
    
    # NK - Green spectrum
    "NK" = "#228B22",
    "NK_CD56bright" = "#90EE90",
    "NK_CD56dim" = "#006400",
    "NK_proliferating" = "#32CD32",
    "NK_adaptive" = "#2E8B57",
    
    # Other T cells - Purple/Brown
    "NKT" = "#9932CC",
    "MAIT" = "#FF8C00",
    "gdT" = "#8B4513",
    "gdT_Vd2" = "#A0522D",
    "gdT_nonVd2" = "#D2691E",
    "dnT" = "#BC8F8F",
    
    # Generic states
    "T_naive" = "#FFB6C1",
    "T_cm" = "#F08080",
    "T_em" = "#CD5C5C",
    "T_proliferating" = "#00CED1",
    "T_activated" = "#FF69B4",
    "T_tissue_resident" = "#DEB887",
    "ISG_high" = "#9400D3",
    
    # Other/Unknown
    "Other" = "#A9A9A9",
    "Unknown" = "#D3D3D3",
    "Unassigned" = "#C0C0C0",
    "Non_T" = "#808080"
  )
}

#' Get paired palette for categorical variables
get_paired_palette <- function(n) {
  if (n <= 12) {
    return(brewer.pal(max(3, n), "Paired")[1:n])
  } else {
    return(colorRampPalette(brewer.pal(12, "Paired"))(n))
  }
}

# ==============================================================================
# 1. MARKER LOADING AND LABEL HARMONIZATION
# ==============================================================================

#' Load T cell marker definitions from YAML
load_tcell_markers <- function(yaml_path = "tcell_markers.yaml") {
  if (!file.exists(yaml_path)) {
    stop("Marker definition file not found: ", yaml_path)
  }
  yaml::read_yaml(yaml_path)
}

#' Get all unique markers from config
get_all_markers <- function(markers_config) {
  all_markers <- c()
  categories <- c("lineages", "cd4_subtypes", "cd8_subtypes", 
                  "other_states", "special", "nk_subtypes")
  
  for (cat in categories) {
    if (cat %in% names(markers_config)) {
      for (cell_type in names(markers_config[[cat]])) {
        ct_info <- markers_config[[cat]][[cell_type]]
        all_markers <- c(all_markers, ct_info$positive, ct_info$negative, ct_info$optional)
      }
    }
  }
  unique(all_markers[!is.na(all_markers)])
}

#' Comprehensive label harmonization map
#' Maps ALL known labels from all sources to canonical names
get_label_harmonization_map <- function() {
  list(
    # =========================================================================
    # Azimuth predicted.celltype.l1
    # =========================================================================
    "CD4 T" = "CD4_T",
    "CD8 T" = "CD8_T",
    "other T" = "T_other",
    "NK" = "NK",
    "B" = "Non_T",
    "Mono" = "Non_T",
    "DC" = "Non_T",
    "other" = "Other",
    
    # =========================================================================
    # Azimuth predicted.celltype.l2
    # =========================================================================
    "CD4 Naive" = "CD4_naive",
    "CD4 TCM" = "CD4_Tcm",
    "CD4 TEM" = "CD4_Tem",
    "CD4 CTL" = "CD4_CTL",
    "CD4 Proliferating" = "CD4_proliferating",
    "CD8 Naive" = "CD8_naive",
    "CD8 TCM" = "CD8_Tcm",
    "CD8 TEM" = "CD8_Tem",
    "CD8 Proliferating" = "CD8_proliferating",
    "Treg" = "Treg",
    "MAIT" = "MAIT",
    "gdT" = "gdT",
    "dnT" = "dnT",
    "NK_CD56bright" = "NK_CD56bright",
    "NK Proliferating" = "NK_proliferating",
    "ILC" = "Other",
    
    # =========================================================================
    # Azimuth predicted.celltype.l3
    # =========================================================================
    "CD4 Naive_1" = "CD4_naive",
    "CD4 Naive_2" = "CD4_naive",
    "CD4 TCM_1" = "CD4_Tcm",
    "CD4 TCM_2" = "CD4_Tcm",
    "CD4 TCM_3" = "CD4_Tcm",
    "CD4 TEM_1" = "CD4_Tem",
    "CD4 TEM_2" = "CD4_Tem",
    "CD4 TEM_3" = "CD4_Tem",
    "CD4 TEM_4" = "CD4_Tem",
    "CD8 Naive_1" = "CD8_naive",
    "CD8 Naive_2" = "CD8_naive",
    "CD8 TCM_1" = "CD8_Tcm",
    "CD8 TCM_2" = "CD8_Tcm",
    "CD8 TCM_3" = "CD8_Tcm",
    "CD8 TEM_1" = "CD8_Tem",
    "CD8 TEM_2" = "CD8_Tem",
    "CD8 TEM_3" = "CD8_Tem",
    "CD8 TEM_4" = "CD8_Tem",
    "CD8 TEM_5" = "CD8_Tem",
    "CD8 TEM_6" = "CD8_Tem",
    "Treg Memory" = "Treg_memory",
    "Treg Naive" = "Treg_naive",
    "NK_1" = "NK",
    "NK_2" = "NK",
    "NK_3" = "NK",
    "NK_4" = "NK",
    "gdT_1" = "gdT",
    "gdT_2" = "gdT",
    "gdT_3" = "gdT",
    "gdT_4" = "gdT",
    "MAIT_1" = "MAIT",
    "MAIT_2" = "MAIT",
    "dnT_1" = "dnT",
    "dnT_2" = "dnT",
    
    # =========================================================================
    # Monaco.labels
    # =========================================================================
    "Naive CD4 T cells" = "CD4_naive",
    "Central memory CD4 T cells" = "CD4_Tcm",
    "Effector memory CD4 T cells" = "CD4_Tem",
    "Terminal effector CD4 T cells" = "CD4_Temra",
    "Naive CD8 T cells" = "CD8_naive",
    "Central memory CD8 T cells" = "CD8_Tcm",
    "Effector memory CD8 T cells" = "CD8_Tem",
    "Terminal effector CD8 T cells" = "CD8_Temra",
    "T regulatory cells" = "Treg",
    "Th1 cells" = "Th1",
    "Th2 cells" = "Th2",
    "Th17 cells" = "Th17",
    "Th1/Th17 cells" = "Th1_Th17",
    "Follicular helper T cells" = "Tfh",
    "MAIT cells" = "MAIT",
    "Vd2 gd T cells" = "gdT_Vd2",
    "Non-Vd2 gd T cells" = "gdT_nonVd2",
    "Natural killer cells" = "NK",
    # Non-T cells from Monaco
    "Low-density basophils" = "Non_T",
    "Low-density neutrophils" = "Non_T",
    "Classical monocytes" = "Non_T",
    "Intermediate monocytes" = "Non_T",
    "Non classical monocytes" = "Non_T",
    "Myeloid dendritic cells" = "Non_T",
    "Plasmacytoid dendritic cells" = "Non_T",
    "Progenitor cells" = "Non_T",
    "Naive B cells" = "Non_T",
    "Non-switched memory B cells" = "Non_T",
    "Switched memory B cells" = "Non_T",
    "Exhausted B cells" = "Non_T",
    "Plasmablasts" = "Non_T",
    
    # =========================================================================
    # HPCA.labels (Human Primary Cell Atlas)
    # =========================================================================
    "T_cells" = "T_other",
    "T_cell:CD4+" = "CD4_T",
    "T_cell:CD4+_naive" = "CD4_naive",
    "T_cell:CD4+_Naive" = "CD4_naive",
    "T_cell:CD4+_central_memory" = "CD4_Tcm",
    "T_cell:CD4+_Central_memory" = "CD4_Tcm",
    "T_cell:CD4+_effector_memory" = "CD4_Tem",
    "T_cell:CD4+_Effector_memory" = "CD4_Tem",
    "T_cell:CD8+" = "CD8_T",
    "T_cell:CD8+_naive" = "CD8_naive",
    "T_cell:CD8+_Naive" = "CD8_naive",
    "T_cell:CD8+_Central_memory" = "CD8_Tcm",
    "T_cell:CD8+_central_memory" = "CD8_Tcm",
    "T_cell:CD8+_effector_memory" = "CD8_Tem",
    "T_cell:CD8+_Effector_memory" = "CD8_Tem",
    "T_cell:CD8+_effector_memory_RA" = "CD8_Temra",
    "T_cell:Treg" = "Treg",
    "T_cell:Treg:Naive" = "Treg_naive",
    "T_cell:gamma-delta" = "gdT",
    "T_cell:effector" = "T_em",
    "NK_cell" = "NK",
    "NK_cell:CD56hiCD62L+" = "NK_CD56bright",
    "NK_cell:CD56hi_CD62L+" = "NK_CD56bright",
    "NK_cell:CD56loCD16+" = "NK_CD56dim",
    "NK_cell:CD56lo_CD16+" = "NK_CD56dim",
    # Non-T from HPCA
    "B_cell" = "Non_T",
    "B_cell:Naive" = "Non_T",
    "B_cell:Memory" = "Non_T",
    "B_cell:Plasma_cell" = "Non_T",
    "Monocyte" = "Non_T",
    "Monocyte:CD14+" = "Non_T",
    "Monocyte:CD16+" = "Non_T",
    "DC" = "Non_T",
    "DC:monocyte-derived" = "Non_T",
    "Macrophage" = "Non_T",
    "Neutrophil" = "Non_T",
    "Erythrocyte" = "Non_T",
    "Platelet" = "Non_T",
    "HSC" = "Non_T",
    "CMP" = "Non_T",
    "GMP" = "Non_T",
    "MEP" = "Non_T",
    
    # =========================================================================
    # Gate subtypes (from marker-based gating) - already canonical
    # =========================================================================
    "CD4_unassigned" = "CD4_T",
    "CD8_unassigned" = "CD8_T",
    "NK_unassigned" = "NK",
    "Unassigned" = "Unknown",
    "Unknown" = "Unknown"
  )
}

#' Harmonize a single label to canonical name
harmonize_label <- function(label, harmonization_map = NULL) {
  if (is.null(harmonization_map)) {
    harmonization_map <- get_label_harmonization_map()
  }
  
  if (is.na(label) || label == "" || is.null(label)) {
    return(NA_character_)
  }
  
  # Direct match
  if (label %in% names(harmonization_map)) {
    return(harmonization_map[[label]])
  }
  
  # Try trimmed version
  label_trimmed <- trimws(label)
  if (label_trimmed %in% names(harmonization_map)) {
    return(harmonization_map[[label_trimmed]])
  }
  
  # Return original if no match found (it might already be canonical)
  label
}

#' Harmonize a vector of labels
harmonize_labels <- function(labels, harmonization_map = NULL) {
  if (is.null(harmonization_map)) {
    harmonization_map <- get_label_harmonization_map()
  }
  sapply(labels, function(l) harmonize_label(l, harmonization_map), USE.NAMES = FALSE)
}

# ==============================================================================
# 2. MARKER EXPRESSION SCORING
# ==============================================================================

#' Calculate marker module scores for each cell
calculate_marker_scores <- function(obj, markers_config, assay = "RNA") {
  
  log_message("Calculating marker module scores...")
  DefaultAssay(obj) <- assay
  
  all_markers <- get_all_markers(markers_config)
  available_markers <- intersect(all_markers, rownames(obj))
  
  log_message("  Available markers: ", length(available_markers), " of ", length(all_markers))
  
  # Extract expression matrix
  expr_mat <- NULL
  tryCatch({
    expr_mat <- LayerData(obj, assay = assay, layer = "data")
    log_message("  Using LayerData for expression extraction")
  }, error = function(e) {
    tryCatch({
      expr_mat <<- GetAssayData(obj, assay = assay, layer = "data")
      log_message("  Using GetAssayData for expression extraction")
    }, error = function(e2) {
      expr_mat <<- GetAssayData(obj, assay = assay, slot = "data")
      log_message("  Using GetAssayData (slot) for expression extraction")
    })
  })
  
  if (is.null(expr_mat)) stop("Could not extract expression matrix")
  
  log_message("  Expression matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " cells")
  
  categories <- c("lineages", "cd4_subtypes", "cd8_subtypes", "other_states", "nk_subtypes")
  scores_calculated <- 0
  
  for (cat in categories) {
    if (!cat %in% names(markers_config)) next
    
    for (cell_type in names(markers_config[[cat]])) {
      ct_info <- markers_config[[cat]][[cell_type]]
      pos_markers <- intersect(ct_info$positive, available_markers)
      score_name <- paste0("Score_", cell_type)
      
      if (length(pos_markers) == 0) {
        obj@meta.data[[score_name]] <- 0
        next
      }
      
      tryCatch({
        marker_idx <- which(rownames(expr_mat) %in% pos_markers)
        marker_expr <- expr_mat[marker_idx, , drop = FALSE]
        
        # Handle S4/sparse matrices
        if (isS4(marker_expr)) marker_expr <- as.matrix(marker_expr)
        if (!is.matrix(marker_expr)) marker_expr <- as.matrix(marker_expr)
        mode(marker_expr) <- "numeric"
        
        # Scale each marker and average
        scaled_expr <- apply(marker_expr, 1, function(x) {
          x <- as.numeric(x)
          rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
          if (rng > 0) (x - min(x, na.rm = TRUE)) / rng else rep(0, length(x))
        })
        
        if (is.matrix(scaled_expr)) {
          obj@meta.data[[score_name]] <- rowMeans(scaled_expr, na.rm = TRUE)
        } else {
          obj@meta.data[[score_name]] <- as.numeric(scaled_expr)
        }
        scores_calculated <- scores_calculated + 1
      }, error = function(e) {
        obj@meta.data[[score_name]] <<- 0
        log_message("    Warning: Scoring failed for ", cell_type)
      })
    }
  }
  
  log_message("  Calculated ", scores_calculated, " cell type scores")
  obj
}

# ==============================================================================
# 3. HIERARCHICAL MARKER-BASED GATING
# ==============================================================================

#' Apply hierarchical gating to cells
apply_hierarchical_gating <- function(obj, markers_config, assay = "RNA") {
  
  log_message("Applying hierarchical marker-based gating...")
  DefaultAssay(obj) <- assay
  
  scoring_params <- markers_config$scoring
  if (is.null(scoring_params)) {
    scoring_params <- list(
      positive_weight = 1.0,
      negative_weight = 0.5,
      min_positive_fraction = 0.4
    )
  }
  
  all_markers <- get_all_markers(markers_config)
  available_markers <- intersect(all_markers, rownames(obj))
  
  # Extract expression for markers only
  expr_mat <- NULL
  tryCatch({
    expr_mat <- LayerData(obj, assay = assay, layer = "data")
  }, error = function(e) {
    tryCatch({
      expr_mat <<- GetAssayData(obj, assay = assay, layer = "data")
    }, error = function(e2) {
      expr_mat <<- GetAssayData(obj, assay = assay, slot = "data")
    })
  })
  
  marker_idx <- which(rownames(expr_mat) %in% available_markers)
  expr_mat <- expr_mat[marker_idx, , drop = FALSE]
  
  # Convert S4 to regular matrix
  if (isS4(expr_mat)) {
    log_message("  Converting to dense matrix...")
    expr_mat <- as.matrix(expr_mat)
  }
  mode(expr_mat) <- "numeric"
  
  available_markers <- rownames(expr_mat)
  n_cells <- ncol(expr_mat)
  
  log_message("  Expression matrix: ", nrow(expr_mat), " markers x ", n_cells, " cells")
  
  # Calculate expression thresholds (median of non-zero values)
  log_message("  Calculating expression thresholds...")
  expr_percentiles <- sapply(available_markers, function(g) {
    x <- expr_mat[g, ]
    x_pos <- x[x > 0]
    if (length(x_pos) > 100) {
      quantile(x_pos, 0.5, na.rm = TRUE)
    } else if (length(x_pos) > 0) {
      quantile(x_pos, 0.25, na.rm = TRUE)
    } else {
      0
    }
  })
  expr_percentiles[is.na(expr_percentiles)] <- 0
  
  # Vectorized scoring function
  score_cell_type <- function(ct_info) {
    pos_markers <- intersect(ct_info$positive, available_markers)
    neg_markers <- intersect(ct_info$negative, available_markers)
    
    if (length(pos_markers) == 0) return(rep(0, n_cells))
    
    # Positive score
    pos_expr <- sapply(pos_markers, function(m) {
      thresh <- expr_percentiles[m]
      if (is.na(thresh)) thresh <- 0
      as.numeric(expr_mat[m, ] > thresh)
    })
    if (is.vector(pos_expr)) pos_expr <- matrix(pos_expr, ncol = 1)
    pos_expr[is.na(pos_expr)] <- 0
    pos_score <- rowMeans(pos_expr, na.rm = TRUE)
    pos_score[is.na(pos_score)] <- 0
    
    # Negative score
    neg_score <- rep(0, n_cells)
    if (length(neg_markers) > 0) {
      neg_expr <- sapply(neg_markers, function(m) {
        thresh <- expr_percentiles[m]
        if (is.na(thresh)) thresh <- 0
        as.numeric(expr_mat[m, ] > thresh)
      })
      if (is.vector(neg_expr)) neg_expr <- matrix(neg_expr, ncol = 1)
      neg_expr[is.na(neg_expr)] <- 0
      neg_score <- rowMeans(neg_expr, na.rm = TRUE)
      neg_score[is.na(neg_score)] <- 0
    }
    
    # Combined score with penalty
    final_score <- pos_score * scoring_params$positive_weight - 
      neg_score * scoring_params$negative_weight
    
    # Penalize low positive fraction
    low_pos <- pos_score < scoring_params$min_positive_fraction
    final_score[low_pos] <- final_score[low_pos] * 0.5
    final_score[is.na(final_score)] <- 0
    
    final_score
  }
  
  # Initialize columns
  obj$gate_lineage <- NA_character_
  obj$gate_lineage_score <- 0
  obj$gate_subtype <- NA_character_
  obj$gate_subtype_score <- 0
  
  # --- Level 1: Lineage Assignment ---
  log_message("  Level 1: Assigning lineages (", n_cells, " cells)...")
  lineages <- markers_config$lineages
  lineage_names <- names(lineages)
  
  lineage_scores <- sapply(lineage_names, function(lin) score_cell_type(lineages[[lin]]))
  if (is.vector(lineage_scores)) {
    lineage_scores <- matrix(lineage_scores, nrow = n_cells)
  }
  colnames(lineage_scores) <- lineage_names
  
  obj$gate_lineage <- apply(lineage_scores, 1, function(row) {
    row[is.na(row)] <- -Inf
    max_val <- max(row)
    if (is.finite(max_val) && max_val > 0) lineage_names[which.max(row)] else "Unknown"
  })
  obj$gate_lineage_score <- apply(lineage_scores, 1, function(row) {
    row[is.na(row)] <- 0
    max(row)
  })
  
  log_message("  Lineage distribution:")
  print(table(obj$gate_lineage))
  
  # --- Level 2: Subtype Assignment ---
  log_message("  Level 2: Assigning subtypes...")
  
  # CD4 subtypes
  cd4_cells <- which(obj$gate_lineage == "CD4_T")
  if (length(cd4_cells) > 0 && "cd4_subtypes" %in% names(markers_config)) {
    log_message("    CD4 subtypes for ", length(cd4_cells), " cells...")
    cd4_subtypes <- markers_config$cd4_subtypes
    cd4_names <- names(cd4_subtypes)
    
    subtype_scores <- sapply(cd4_names, function(st) {
      score_cell_type(cd4_subtypes[[st]])[cd4_cells]
    })
    if (is.vector(subtype_scores)) {
      subtype_scores <- matrix(subtype_scores, nrow = length(cd4_cells))
    }
    colnames(subtype_scores) <- cd4_names
    
    assignments <- apply(subtype_scores, 1, function(row) {
      row[is.na(row)] <- -Inf
      max_val <- max(row)
      if (is.finite(max_val) && max_val > 0.1) cd4_names[which.max(row)] else "CD4_T"
    })
    scores <- apply(subtype_scores, 1, function(row) {
      row[is.na(row)] <- 0
      max(row)
    })
    
    obj$gate_subtype[cd4_cells] <- assignments
    obj$gate_subtype_score[cd4_cells] <- scores
  }
  
  # CD8 subtypes
  cd8_cells <- which(obj$gate_lineage == "CD8_T")
  if (length(cd8_cells) > 0 && "cd8_subtypes" %in% names(markers_config)) {
    log_message("    CD8 subtypes for ", length(cd8_cells), " cells...")
    cd8_subtypes <- markers_config$cd8_subtypes
    cd8_names <- names(cd8_subtypes)
    
    subtype_scores <- sapply(cd8_names, function(st) {
      score_cell_type(cd8_subtypes[[st]])[cd8_cells]
    })
    if (is.vector(subtype_scores)) {
      subtype_scores <- matrix(subtype_scores, nrow = length(cd8_cells))
    }
    colnames(subtype_scores) <- cd8_names
    
    assignments <- apply(subtype_scores, 1, function(row) {
      row[is.na(row)] <- -Inf
      max_val <- max(row)
      if (is.finite(max_val) && max_val > 0.1) cd8_names[which.max(row)] else "CD8_T"
    })
    scores <- apply(subtype_scores, 1, function(row) {
      row[is.na(row)] <- 0
      max(row)
    })
    
    obj$gate_subtype[cd8_cells] <- assignments
    obj$gate_subtype_score[cd8_cells] <- scores
  }
  
  # NK subtypes
  nk_cells <- which(obj$gate_lineage == "NK")
  if (length(nk_cells) > 0 && "nk_subtypes" %in% names(markers_config)) {
    log_message("    NK subtypes for ", length(nk_cells), " cells...")
    nk_subtypes <- markers_config$nk_subtypes
    nk_names <- names(nk_subtypes)
    
    subtype_scores <- sapply(nk_names, function(st) {
      score_cell_type(nk_subtypes[[st]])[nk_cells]
    })
    if (is.vector(subtype_scores)) {
      subtype_scores <- matrix(subtype_scores, nrow = length(nk_cells))
    }
    colnames(subtype_scores) <- nk_names
    
    assignments <- apply(subtype_scores, 1, function(row) {
      row[is.na(row)] <- -Inf
      max_val <- max(row)
      if (is.finite(max_val) && max_val > 0.1) nk_names[which.max(row)] else "NK"
    })
    scores <- apply(subtype_scores, 1, function(row) {
      row[is.na(row)] <- 0
      max(row)
    })
    
    obj$gate_subtype[nk_cells] <- assignments
    obj$gate_subtype_score[nk_cells] <- scores
  }
  
  # Other lineages keep their lineage as subtype
  other_lineages <- c("gdT", "MAIT", "NKT", "dnT")
  for (lin in other_lineages) {
    lin_cells <- which(obj$gate_lineage == lin)
    if (length(lin_cells) > 0) {
      obj$gate_subtype[lin_cells] <- lin
      obj$gate_subtype_score[lin_cells] <- obj$gate_lineage_score[lin_cells]
    }
  }
  
  # Fill remaining NAs
  obj$gate_subtype[is.na(obj$gate_subtype)] <- "Unknown"
  
  log_message("  Gating complete")
  obj
}

# ==============================================================================
# 4. COMPREHENSIVE CONSENSUS ANNOTATION
# ==============================================================================

#' Get majority vote per cluster with harmonized labels
get_cluster_majority_harmonized <- function(obj, annotation_col, cluster_col,
                                            harmonization_map, min_fraction = 0.15) {
  
  if (!annotation_col %in% colnames(obj@meta.data)) {
    return(NULL)
  }
  
  # Get and harmonize labels
  labels <- obj@meta.data[[annotation_col]]
  harmonized <- harmonize_labels(labels, harmonization_map)
  
  df <- data.frame(
    cluster = obj@meta.data[[cluster_col]],
    label = harmonized,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(label) & label != "" & label != "Unknown" & label != "Non_T") %>%
    group_by(cluster, label) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(
      total = sum(count),
      fraction = count / total
    ) %>%
    slice_max(count, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  df$confident <- df$fraction >= min_fraction
  df
}

#' Generate comprehensive cluster consensus using ALL annotation sources
generate_cluster_consensus <- function(obj, 
                                       cluster_col = "leiden.cluster",
                                       markers_config = NULL) {
  
  log_message("Generating comprehensive cluster consensus...")
  
  harmonization_map <- get_label_harmonization_map()
  
  # Define ALL annotation sources with EQUAL weights
  # Each source contributes equally to the consensus
  annotation_sources <- list(
    # Marker-based gating (expression modules)
    list(col = "gate_subtype", weight = 1.0, name = "Gate_Subtype"),
    
    # Azimuth predictions (all levels equal)
    list(col = "predicted.celltype.l3", weight = 1.0, name = "Azimuth_L3"),
    list(col = "predicted.celltype.l2", weight = 1.0, name = "Azimuth_L2"),
    list(col = "predicted.celltype.l1", weight = 1.0, name = "Azimuth_L1"),
    
    # SingleR reference datasets (equal weight)
    list(col = "Monaco.labels", weight = 1.0, name = "Monaco"),
    list(col = "HPCA.labels", weight = 1.0, name = "HPCA")
  )
  
  clusters <- sort(unique(obj@meta.data[[cluster_col]]))
  log_message("  Found ", length(clusters), " clusters")
  
  # Initialize consensus dataframe
  consensus_df <- data.frame(cluster = clusters, stringsAsFactors = FALSE)
  
  # Track which sources are available
  sources_used <- c()
  total_weight <- 0
  
  # Process each annotation source
  for (source in annotation_sources) {
    col_name <- source$col
    weight <- source$weight
    display_name <- source$name
    
    if (!col_name %in% colnames(obj@meta.data)) {
      log_message("    Skipping: ", col_name, " (not found)")
      next
    }
    
    # Check if column has valid data
    valid_data <- sum(!is.na(obj@meta.data[[col_name]])) > 0
    if (!valid_data) {
      log_message("    Skipping: ", col_name, " (no valid data)")
      next
    }
    
    log_message("    Processing: ", col_name, " (weight: ", weight, ")")
    sources_used <- c(sources_used, display_name)
    total_weight <- total_weight + weight
    
    maj <- get_cluster_majority_harmonized(obj, col_name, cluster_col, harmonization_map)
    
    if (is.null(maj) || nrow(maj) == 0) {
      consensus_df[[paste0(display_name, "_label")]] <- NA_character_
      consensus_df[[paste0(display_name, "_fraction")]] <- NA_real_
    } else {
      match_idx <- match(as.character(consensus_df$cluster), as.character(maj$cluster))
      consensus_df[[paste0(display_name, "_label")]] <- maj$label[match_idx]
      consensus_df[[paste0(display_name, "_fraction")]] <- maj$fraction[match_idx]
    }
    
    consensus_df[[paste0(display_name, "_weight")]] <- weight
  }
  
  log_message("    Sources used: ", paste(sources_used, collapse = ", "))
  log_message("    Total weight: ", round(total_weight, 2))
  
  # Calculate cluster sizes
  cluster_sizes <- table(obj@meta.data[[cluster_col]])
  consensus_df$n_cells <- as.numeric(cluster_sizes[as.character(consensus_df$cluster)])
  
  # Assign final annotations using weighted voting
  consensus_df <- assign_final_annotation_weighted(consensus_df, sources_used)
  
  consensus_df
}

#' Assign final annotation using weighted voting across all sources
assign_final_annotation_weighted <- function(consensus_df, sources_used) {
  
  log_message("  Assigning final annotations with weighted voting...")
  
  # Get label and weight columns
  label_cols <- paste0(sources_used, "_label")
  fraction_cols <- paste0(sources_used, "_fraction")
  weight_cols <- paste0(sources_used, "_weight")
  
  final_annotations <- character(nrow(consensus_df))
  annotation_confidence <- numeric(nrow(consensus_df))
  n_sources_agree <- integer(nrow(consensus_df))
  
  for (i in 1:nrow(consensus_df)) {
    row <- consensus_df[i, ]
    
    # Collect weighted votes for each unique label
    votes <- list()
    
    for (j in seq_along(sources_used)) {
      src <- sources_used[j]
      label_col <- label_cols[j]
      frac_col <- fraction_cols[j]
      weight_col <- weight_cols[j]
      
      if (!label_col %in% names(row)) next
      
      label <- row[[label_col]]
      if (is.na(label) || label == "" || label == "Unknown" || label == "Non_T") next
      
      frac <- ifelse(!is.na(row[[frac_col]]), row[[frac_col]], 0.5)
      weight <- ifelse(!is.na(row[[weight_col]]), row[[weight_col]], 0.1)
      
      # Weighted score = source_weight * within_cluster_fraction
      weighted_score <- weight * frac
      
      if (label %in% names(votes)) {
        votes[[label]] <- votes[[label]] + weighted_score
      } else {
        votes[[label]] <- weighted_score
      }
    }
    
    if (length(votes) == 0) {
      final_annotations[i] <- "Unknown"
      annotation_confidence[i] <- 0
      n_sources_agree[i] <- 0
    } else {
      # Best label = highest weighted score
      best_label <- names(votes)[which.max(unlist(votes))]
      final_annotations[i] <- best_label
      
      # Confidence = sum of weights from agreeing sources / total weight used
      agreeing_weight <- 0
      n_agree <- 0
      for (j in seq_along(sources_used)) {
        label_col <- label_cols[j]
        frac_col <- fraction_cols[j]
        weight_col <- weight_cols[j]
        
        if (!label_col %in% names(row)) next
        label <- row[[label_col]]
        if (is.na(label)) next
        
        if (label == best_label) {
          frac <- ifelse(!is.na(row[[frac_col]]), row[[frac_col]], 0.5)
          weight <- ifelse(!is.na(row[[weight_col]]), row[[weight_col]], 0.1)
          agreeing_weight <- agreeing_weight + weight * frac
          n_agree <- n_agree + 1
        }
      }
      
      annotation_confidence[i] <- agreeing_weight
      n_sources_agree[i] <- n_agree
    }
  }
  
  consensus_df$final_annotation <- final_annotations
  consensus_df$annotation_confidence <- annotation_confidence
  consensus_df$n_sources_agree <- n_sources_agree
  
  # Report confidence distribution
  log_message("  Confidence distribution:")
  log_message("    Low (<0.3): ", sum(annotation_confidence < 0.3), " clusters")
  log_message("    Medium (0.3-0.5): ", sum(annotation_confidence >= 0.3 & annotation_confidence < 0.5), " clusters")
  log_message("    High (>=0.5): ", sum(annotation_confidence >= 0.5), " clusters")
  
  consensus_df
}

#' Apply cluster annotations to Seurat object
apply_cluster_annotations <- function(obj, consensus_df, cluster_col = "leiden.cluster") {
  
  log_message("Applying cluster annotations to cells...")
  
  cell_clusters <- as.character(obj@meta.data[[cluster_col]])
  
  # Final annotation
  annotation_map <- setNames(consensus_df$final_annotation, 
                             as.character(consensus_df$cluster))
  obj@meta.data$cluster_annotation <- unname(annotation_map[cell_clusters])
  obj@meta.data$cluster_annotation[is.na(obj@meta.data$cluster_annotation)] <- "Unknown"
  
  # Confidence
  confidence_map <- setNames(consensus_df$annotation_confidence,
                             as.character(consensus_df$cluster))
  obj@meta.data$annotation_confidence <- unname(confidence_map[cell_clusters])
  obj@meta.data$annotation_confidence[is.na(obj@meta.data$annotation_confidence)] <- 0
  
  # N sources agree
  agree_map <- setNames(consensus_df$n_sources_agree,
                        as.character(consensus_df$cluster))
  obj@meta.data$n_sources_agree <- unname(agree_map[cell_clusters])
  obj@meta.data$n_sources_agree[is.na(obj@meta.data$n_sources_agree)] <- 0
  
  log_message("  Annotations applied to ", ncol(obj), " cells")
  log_message("  Annotation distribution:")
  print(table(obj@meta.data$cluster_annotation))
  
  obj
}

# ==============================================================================
# 5. PER-CLUSTER STATISTICS
# ==============================================================================

#' Get key canonical markers for analysis
get_key_markers <- function() {
  c(
    # Lineage markers
    "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "NCAM1",
    # Naive/Memory
    "CCR7", "SELL", "LEF1", "TCF7", "IL7R", "CD27", "CD28",
    # Effector/Cytotoxic
    "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "FGFBP2", "KLRG1",
    # CD4 Th subsets
    "TBX21", "IFNG", "GATA3", "RORC", "IL17A", "BCL6", "CXCR5",
    # Treg
    "FOXP3", "IL2RA", "CTLA4", "IKZF2",
    # Exhaustion/Activation
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CD69", "CD38",
    # Other
    "KLRB1", "TRDC", "SLC4A10", "MKI67", "TOP2A"
  )
}

#' Calculate comprehensive per-cluster marker statistics
calculate_cluster_marker_stats <- function(obj, markers = NULL, 
                                           cluster_col = "cluster_annotation",
                                           assay = "RNA") {
  
  log_message("Calculating per-cluster marker statistics...")
  
  if (is.null(markers)) {
    markers <- get_key_markers()
  }
  
  DefaultAssay(obj) <- assay
  available_markers <- intersect(markers, rownames(obj))
  
  if (length(available_markers) == 0) {
    warning("No markers found in object")
    return(NULL)
  }
  
  log_message("  Using ", length(available_markers), " markers")
  
  # Extract expression
  expr_mat <- NULL
  tryCatch({
    expr_mat <- LayerData(obj, assay = assay, layer = "data")
  }, error = function(e) {
    tryCatch({
      expr_mat <<- GetAssayData(obj, assay = assay, layer = "data")
    }, error = function(e2) {
      expr_mat <<- GetAssayData(obj, assay = assay, slot = "data")
    })
  })
  
  marker_idx <- which(rownames(expr_mat) %in% available_markers)
  expr_mat <- expr_mat[marker_idx, , drop = FALSE]
  
  if (isS4(expr_mat)) expr_mat <- as.matrix(expr_mat)
  mode(expr_mat) <- "numeric"
  
  clusters <- obj@meta.data[[cluster_col]]
  cluster_levels <- sort(unique(clusters))
  
  # Calculate statistics for each cluster
  stats_list <- list()
  
  for (cl in cluster_levels) {
    cells <- which(clusters == cl)
    if (length(cells) == 0) next
    
    cl_expr <- expr_mat[, cells, drop = FALSE]
    
    for (marker in rownames(cl_expr)) {
      vals <- cl_expr[marker, ]
      
      stats_list[[paste(cl, marker, sep = "_")]] <- data.frame(
        cluster = cl,
        marker = marker,
        mean_expression = mean(vals, na.rm = TRUE),
        median_expression = median(vals, na.rm = TRUE),
        pct_expressing = mean(vals > 0, na.rm = TRUE) * 100,
        sd_expression = sd(vals, na.rm = TRUE),
        n_cells = length(cells),
        stringsAsFactors = FALSE
      )
    }
  }
  
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL
  
  log_message("  Calculated stats for ", length(available_markers), 
              " markers across ", length(cluster_levels), " clusters")
  
  stats_df
}

#' Generate per-cluster summary with key marker expression
generate_cluster_summary <- function(obj, consensus_df, 
                                     cluster_col = "leiden.cluster",
                                     key_markers = NULL) {
  
  log_message("Generating cluster summary with marker expression...")
  
  if (is.null(key_markers)) {
    # Use subset of most informative markers
    key_markers <- c("CD4", "CD8A", "FOXP3", "CCR7", "GZMB", "PDCD1", 
                     "NCAM1", "KLRB1", "MKI67")
  }
  
  available <- intersect(key_markers, rownames(obj))
  
  if (length(available) > 0) {
    # Calculate mean expression per cluster
    stats <- calculate_cluster_marker_stats(obj, markers = available, 
                                            cluster_col = cluster_col)
    
    # Pivot to wide format
    if (!is.null(stats)) {
      mean_wide <- stats %>%
        select(cluster, marker, mean_expression) %>%
        pivot_wider(names_from = marker, values_from = mean_expression,
                    names_prefix = "mean_")
      
      pct_wide <- stats %>%
        select(cluster, marker, pct_expressing) %>%
        pivot_wider(names_from = marker, values_from = pct_expressing,
                    names_prefix = "pct_")
      
      # Merge with consensus
      summary_df <- consensus_df %>%
        left_join(mean_wide, by = "cluster") %>%
        left_join(pct_wide, by = "cluster")
      
      return(summary_df)
    }
  }
  
  consensus_df
}

# ==============================================================================
# 6. PUBLICATION-QUALITY VISUALIZATIONS
# ==============================================================================

#' Create UMAP plots comparing annotations
plot_annotation_umaps <- function(obj, reduction = "umap", pt.size = 0.3,
                                  output_dir = NULL) {
  
  log_message("Creating annotation UMAP plots...")
  
  # Columns to plot
  plot_cols <- c("cluster_annotation", "leiden.cluster", "gate_lineage",
                 "predicted.celltype.l2", "Monaco.labels", "HPCA.labels")
  plot_cols <- plot_cols[plot_cols %in% colnames(obj@meta.data)]
  
  plots <- list()
  
  for (col in plot_cols) {
    # Clean title
    title <- gsub("\\.", " ", col)
    title <- gsub("_", " ", title)
    title <- tools::toTitleCase(title)
    
    vals <- unique(na.omit(obj@meta.data[[col]]))
    n_vals <- length(vals)
    
    # Use Paired palette for all categorical variables
    if (n_vals <= 12) {
      cols <- brewer.pal(max(3, n_vals), "Paired")[1:n_vals]
    } else {
      cols <- colorRampPalette(brewer.pal(12, "Paired"))(n_vals)
    }
    names(cols) <- vals
    
    p <- DimPlot(obj, reduction = reduction, group.by = col,
                 cols = cols, pt.size = pt.size,
                 label = TRUE, label.size = 3, repel = TRUE) +
      ggtitle(title) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    plots[[col]] <- p
  }
  
  # Combine
  n_plots <- length(plots)
  ncol <- min(3, n_plots)
  nrow <- ceiling(n_plots / ncol)
  
  combined <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = "Cell Type Annotation Comparison",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, "08_annotation_umaps.pdf"), combined,
           width = 6 * ncol, height = 6 * nrow, dpi = 300)
    ggsave(file.path(output_dir, "08_annotation_umaps.png"), combined,
           width = 6 * ncol, height = 6 * nrow, dpi = 300)
    log_message("  Saved: annotation_umaps.pdf/png")
  }
  
  combined
}

#' Create canonical marker DotPlot with viridis palette
plot_canonical_dotplot <- function(obj, group.by = "cluster_annotation",
                                   output_dir = NULL) {
  
  log_message("Creating canonical marker DotPlot...")
  
  key_markers <- get_key_markers()
  available <- intersect(key_markers, rownames(obj))
  
  if (length(available) == 0) {
    warning("No key markers found")
    return(NULL)
  }
  
  # Order markers by category
  marker_order <- c(
    "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "NCAM1",  # Lineage
    "CCR7", "SELL", "LEF1", "TCF7", "IL7R",          # Naive
    "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7",  # Cytotoxic
    "TBX21", "IFNG", "GATA3", "RORC", "BCL6",        # Th TFs
    "FOXP3", "IL2RA", "CTLA4",                        # Treg
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX",       # Exhaustion
    "KLRB1", "TRDC", "MKI67"                          # Other
  )
  available_ordered <- marker_order[marker_order %in% available]
  available_ordered <- c(available_ordered, setdiff(available, available_ordered))
  
  p <- DotPlot(obj, features = available_ordered, group.by = group.by) +
    scale_color_viridis_c(option = "viridis", name = "Avg Expr") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(title = "Canonical T Cell Marker Expression",
         size = "% Expr")
  
  if (!is.null(output_dir)) {
    n_groups <- length(unique(obj@meta.data[[group.by]]))
    height <- max(8, n_groups * 0.4)
    width <- max(12, length(available_ordered) * 0.35)
    
    ggsave(file.path(output_dir, "08_marker_dotplot.pdf"), p,
           width = width, height = height, dpi = 300)
    ggsave(file.path(output_dir, "08_marker_dotplot.png"), p,
           width = width, height = height, dpi = 300)
    log_message("  Saved: marker_dotplot.pdf/png")
  }
  
  p
}

#' Create key marker FeaturePlots with viridis palette
plot_key_features <- function(obj, reduction = "umap", output_dir = NULL) {
  
  log_message("Creating key marker FeaturePlots...")
  
  markers <- c("CD4", "CD8A", "FOXP3", "GZMB", "CCR7", "PDCD1",
               "NCAM1", "KLRB1", "MKI67")
  available <- intersect(markers, rownames(obj))
  
  if (length(available) == 0) {
    warning("No key markers found")
    return(NULL)
  }
  
  # Create individual plots with better visibility
  plots <- list()
  for (marker in available) {
    p <- FeaturePlot(obj, features = marker, reduction = reduction,
                     order = TRUE, pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q95") +
      scale_color_viridis_c(option = "viridis", na.value = "grey85") +
      ggtitle(marker) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.8, "cm"),
        legend.title = element_text(size = 9),
        plot.margin = margin(5, 5, 5, 5)
      )
    plots[[marker]] <- p
  }
  
  # Combine plots
  combined <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = "Key T Cell Marker Expression",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  if (!is.null(output_dir)) {
    nrow <- ceiling(length(available) / 3)
    ggsave(file.path(output_dir, "08_key_marker_features.pdf"), combined,
           width = 15, height = 5 * nrow, dpi = 300)
    ggsave(file.path(output_dir, "08_key_marker_features.png"), combined,
           width = 15, height = 5 * nrow, dpi = 300)
    log_message("  Saved: key_marker_features.pdf/png")
  }
  
  combined
}

#' Create cluster composition stacked bar plot
plot_cluster_composition <- function(obj, cluster_col = "leiden.cluster",
                                     annotation_col = "cluster_annotation",
                                     output_dir = NULL) {
  
  log_message("Creating cluster composition plot...")
  
  comp_df <- obj@meta.data %>%
    group_by(.data[[cluster_col]], .data[[annotation_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup()
  
  # Use Paired palette for categorical annotations
  annot_levels <- unique(comp_df[[annotation_col]])
  n_annot <- length(annot_levels)
  if (n_annot <= 12) {
    cols <- brewer.pal(max(3, n_annot), "Paired")[1:n_annot]
  } else {
    cols <- colorRampPalette(brewer.pal(12, "Paired"))(n_annot)
  }
  names(cols) <- annot_levels
  
  p <- ggplot(comp_df, aes(x = factor(.data[[cluster_col]]),
                           y = fraction,
                           fill = .data[[annotation_col]])) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    scale_fill_manual(values = cols, name = "Annotation") +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(x = "Cluster", y = "Fraction",
         title = "Cluster Composition by Final Annotation")
  
  if (!is.null(output_dir)) {
    n_clusters <- length(unique(comp_df[[cluster_col]]))
    width <- max(10, n_clusters * 0.6)
    
    ggsave(file.path(output_dir, "08_cluster_composition.pdf"), p,
           width = width, height = 7, dpi = 300)
    ggsave(file.path(output_dir, "08_cluster_composition.png"), p,
           width = width, height = 7, dpi = 300)
    log_message("  Saved: cluster_composition.pdf/png")
  }
  
  p
}

#' Create annotation confidence bar plot
plot_annotation_confidence <- function(consensus_df, output_dir = NULL) {
  
  log_message("Creating annotation confidence plot...")
  
  df <- consensus_df %>%
    mutate(cluster = factor(cluster, 
                            levels = cluster[order(annotation_confidence, decreasing = TRUE)])) %>%
    arrange(desc(annotation_confidence))
  
  # Use Paired palette for annotations
  annot_levels <- unique(df$final_annotation)
  n_annot <- length(annot_levels)
  if (n_annot <= 12) {
    cols <- brewer.pal(max(3, n_annot), "Paired")[1:n_annot]
  } else {
    cols <- colorRampPalette(brewer.pal(12, "Paired"))(n_annot)
  }
  names(cols) <- annot_levels
  
  p <- ggplot(df, aes(x = cluster, y = annotation_confidence, fill = final_annotation)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red", alpha = 0.8, linewidth = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", alpha = 0.8, linewidth = 0.7) +
    scale_fill_manual(values = cols, name = "Annotation") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df$annotation_confidence, na.rm = TRUE) * 1.1)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    ) +
    labs(x = "Cluster", y = "Annotation Confidence Score",
         title = "Per-Cluster Annotation Confidence",
         subtitle = "Red dashed: low threshold (0.3), Green dashed: high threshold (0.5)")
  
  if (!is.null(output_dir)) {
    n_clusters <- nrow(df)
    width <- max(10, n_clusters * 0.5)
    
    ggsave(file.path(output_dir, "08_annotation_confidence.pdf"), p,
           width = width, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "08_annotation_confidence.png"), p,
           width = width, height = 6, dpi = 300)
    log_message("  Saved: annotation_confidence.pdf/png")
  }
  
  p
}

#' Create annotation source agreement heatmap
plot_source_agreement <- function(consensus_df, output_dir = NULL) {
  
  log_message("Creating source agreement heatmap...")
  
  # Get label columns
  label_cols <- grep("_label$", names(consensus_df), value = TRUE)
  
  if (length(label_cols) < 2) {
    log_message("  Not enough sources for agreement plot")
    return(NULL)
  }
  
  # Build matrix
  source_names <- sub("_label$", "", label_cols)
  
  # Create tidy data for plotting
  plot_data <- data.frame()
  for (i in 1:nrow(consensus_df)) {
    for (lc in label_cols) {
      source <- sub("_label$", "", lc)
      label <- consensus_df[i, lc]
      final <- consensus_df$final_annotation[i]
      agrees <- !is.na(label) && label == final
      
      plot_data <- rbind(plot_data, data.frame(
        cluster = consensus_df$cluster[i],
        source = source,
        label = ifelse(is.na(label), "NA", 
                       ifelse(nchar(label) > 10, paste0(substr(label, 1, 9), ".."), label)),
        agrees_with_final = agrees,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Plot as tile with better colors
  p <- ggplot(plot_data, aes(x = source, y = factor(cluster), fill = agrees_with_final)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 2.5, color = "black") +
    scale_fill_manual(values = c("TRUE" = "#66C2A5", "FALSE" = "#FC8D62"),
                      name = "Agrees with\nFinal Annotation",
                      labels = c("TRUE" = "Yes", "FALSE" = "No")) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(title = "Annotation Source Agreement with Final Consensus")
  
  if (!is.null(output_dir)) {
    n_clusters <- nrow(consensus_df)
    n_sources <- length(label_cols)
    
    ggsave(file.path(output_dir, "08_source_agreement.pdf"), p,
           width = max(8, n_sources * 1.5), height = max(8, n_clusters * 0.4), dpi = 300)
    ggsave(file.path(output_dir, "08_source_agreement.png"), p,
           width = max(8, n_sources * 1.5), height = max(8, n_clusters * 0.4), dpi = 300)
    log_message("  Saved: source_agreement.pdf/png")
  }
  
  p
}

#' Create marker heatmap
plot_marker_heatmap_ggplot <- function(obj, cluster_col = "cluster_annotation",
                                       markers = NULL, output_dir = NULL) {
  
  log_message("Creating marker expression heatmap...")
  
  if (is.null(markers)) {
    markers <- get_key_markers()
  }
  
  stats <- calculate_cluster_marker_stats(obj, markers = markers, cluster_col = cluster_col)
  
  if (is.null(stats) || nrow(stats) == 0) {
    log_message("  No stats available for heatmap")
    return(NULL)
  }
  
  # Scale expression within each marker
  scaled_stats <- stats %>%
    group_by(marker) %>%
    mutate(scaled_expr = scale(mean_expression)[,1]) %>%
    ungroup()
  
  # Clip values
  scaled_stats$scaled_expr[scaled_stats$scaled_expr > 2] <- 2
  scaled_stats$scaled_expr[scaled_stats$scaled_expr < -2] <- -2
  
  # Use viridis for heatmap
  p <- ggplot(scaled_stats, aes(x = cluster, y = marker, fill = scaled_expr)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_viridis_c(option = "viridis", limits = c(-2, 2), name = "Z-score") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(title = "Marker Expression by Cluster (Z-scaled)")
  
  if (!is.null(output_dir)) {
    n_clusters <- length(unique(stats$cluster))
    n_markers <- length(unique(stats$marker))
    
    ggsave(file.path(output_dir, "08_marker_heatmap.pdf"), p,
           width = max(10, n_clusters * 0.5), height = max(10, n_markers * 0.3), dpi = 300)
    ggsave(file.path(output_dir, "08_marker_heatmap.png"), p,
           width = max(10, n_clusters * 0.5), height = max(10, n_markers * 0.3), dpi = 300)
    log_message("  Saved: marker_heatmap.pdf/png")
  }
  
  p
}

# ==============================================================================
# 7. MAIN ANNOTATION FUNCTION
# ==============================================================================

#' Run complete T cell annotation pipeline
#'
#' @param obj Seurat object
#' @param markers_yaml Path to marker definition YAML file
#' @param cluster_col Column containing cluster IDs
#' @param output_dir Directory to save outputs
#' @param save_object Whether to save the annotated Seurat object
#' @param create_plots Whether to generate visualization plots
#' @return List containing annotated object, consensus dataframe, and statistics
#'
#' @export
run_tcell_annotation <- function(obj,
                                 markers_yaml = "tcell_markers.yaml",
                                 cluster_col = "leiden.cluster",
                                 output_dir = "results/annotations",
                                 save_object = TRUE,
                                 create_plots = TRUE, 
                                 config = NULL) {
  
  plot_dir <- "figs"
  log_message("=== Starting T Cell Annotation Pipeline ===")
  log_message("Cells: ", ncol(obj))
  log_message("Cluster column: ", cluster_col)
  log_message("Output directory: ", output_dir)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message("Created output directory: ", output_dir)
  }
  
  # Load marker definitions
  log_message("Loading marker definitions from: ", markers_yaml)
  markers_config <- load_tcell_markers(markers_yaml)
  log_message("Loaded marker definitions")
  
  # Report available annotation columns
  log_message("Checking available annotation sources...")
  annot_cols <- c("predicted.celltype.l1", "predicted.celltype.l2", 
                  "predicted.celltype.l3", "Monaco.labels", "HPCA.labels")
  for (col in annot_cols) {
    if (col %in% colnames(obj@meta.data)) {
      n_valid <- sum(!is.na(obj@meta.data[[col]]))
      n_unique <- length(unique(na.omit(obj@meta.data[[col]])))
      log_message("  ", col, ": ", n_valid, " cells, ", n_unique, " unique labels")
    } else {
      log_message("  ", col, ": NOT FOUND")
    }
  }
  
  # Step 1: Calculate marker scores
  obj <- calculate_marker_scores(obj, markers_config, assay = "RNA")
  
  # Step 2: Apply hierarchical gating
  obj <- apply_hierarchical_gating(obj, markers_config, assay = "RNA")
  
  # Step 3: Generate cluster consensus
  consensus_df <- generate_cluster_consensus(
    obj,
    cluster_col = cluster_col,
    markers_config = markers_config
  )
  
  # Step 4: Apply annotations to object
  obj <- apply_cluster_annotations(obj, consensus_df, cluster_col = cluster_col)
  
  # Step 5: Calculate cluster statistics
  log_message("Calculating cluster marker statistics...")
  cluster_stats <- calculate_cluster_marker_stats(
    obj,
    markers = get_key_markers(),
    cluster_col = "cluster_annotation",
    assay = "RNA"
  )
  
  # Step 6: Generate summary with marker info
  summary_df <- generate_cluster_summary(
    obj,
    consensus_df,
    cluster_col = cluster_col
  )
  
  # Save tables
  log_message("Saving summary tables...")
  
  write.csv(consensus_df,
            file.path(output_dir, "cluster_consensus_full.csv"),
            row.names = FALSE)
  log_message("  Saved: cluster_consensus_full.csv")
  
  write.csv(summary_df,
            file.path(output_dir, "cluster_summary.csv"),
            row.names = FALSE)
  log_message("  Saved: cluster_summary.csv")
  
  if (!is.null(cluster_stats)) {
    write.csv(cluster_stats,
              file.path(output_dir, "cluster_marker_stats.csv"),
              row.names = FALSE)
    log_message("  Saved: cluster_marker_stats.csv")
  }
  
  # Export per-cell annotations
  cell_annot <- obj@meta.data[, c("cluster_annotation", "annotation_confidence",
                                  "gate_lineage", "gate_subtype"), drop = FALSE]
  cell_annot$cell_barcode <- rownames(cell_annot)
  cell_annot <- cell_annot[, c("cell_barcode", names(cell_annot)[-ncol(cell_annot)])]
  
  write.csv(cell_annot,
            file.path(output_dir, "cell_annotations.csv"),
            row.names = FALSE)
  log_message("  Saved: cell_annotations.csv")
  
  # Generate annotation summary
  annot_summary <- as.data.frame(table(obj@meta.data$cluster_annotation))
  colnames(annot_summary) <- c("annotation", "n_cells")
  annot_summary$pct <- annot_summary$n_cells / ncol(obj) * 100
  annot_summary <- annot_summary[order(-annot_summary$n_cells), ]
  
  write.csv(annot_summary,
            file.path(output_dir, "annotation_summary.csv"),
            row.names = FALSE)
  log_message("  Saved: annotation_summary.csv")
  
  # Step 7: Create visualizations
  if (create_plots) {
    log_message("Creating visualizations...")
    
    plots <- list()
    
    # UMAPs
    tryCatch({
      plots$umaps <- plot_annotation_umaps(obj, output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: UMAP plot failed - ", e$message)
    })
    
    # DotPlot
    tryCatch({
      plots$dotplot <- plot_canonical_dotplot(obj, group.by = "cluster_annotation",
                                              output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: DotPlot failed - ", e$message)
    })
    
    # FeaturePlots
    tryCatch({
      plots$features <- plot_key_features(obj, output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: FeaturePlot failed - ", e$message)
    })
    
    # Cluster composition
    tryCatch({
      plots$composition <- plot_cluster_composition(obj, cluster_col = cluster_col,
                                                    output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: Composition plot failed - ", e$message)
    })
    
    # Confidence plot
    tryCatch({
      plots$confidence <- plot_annotation_confidence(consensus_df, output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: Confidence plot failed - ", e$message)
    })
    
    # Source agreement
    tryCatch({
      plots$agreement <- plot_source_agreement(consensus_df, output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: Agreement plot failed - ", e$message)
    })
    
    # Marker heatmap
    tryCatch({
      plots$heatmap <- plot_marker_heatmap_ggplot(obj, cluster_col = "cluster_annotation",
                                                  output_dir = plot_dir)
    }, error = function(e) {
      log_message("  Warning: Heatmap failed - ", e$message)
    })
  }
  
  # Step 8: Save annotated object
  if (save_object) {
    log_message("Saving annotated Seurat object...")
    obj_path <- file.path(output_dir, "annotated_object.rds")
    saveRDS(obj, obj_path)
    log_message("  Saved: annotated_object.rds")
  }
  
  # Final summary
  log_message("=== Annotation Pipeline Complete ===")
  log_message("Final annotation distribution:")
  print(table(obj@meta.data$cluster_annotation))
  
  # Return results
  list(
    object = obj,
    consensus = consensus_df,
    summary = summary_df,
    stats = cluster_stats,
    annotation_summary = annot_summary
  )
}

#' Quick validation of annotations
#'
#' @param obj Annotated Seurat object
#' @param markers_config Marker configuration
#' @return Validation dataframe
validate_annotations <- function(obj, markers_config) {
  
  log_message("Validating annotations...")
  
  annotations <- unique(obj@meta.data$cluster_annotation)
  annotations <- annotations[!is.na(annotations) & annotations != "Unknown"]
  
  harmonization <- get_label_harmonization_map()
  all_markers <- get_all_markers(markers_config)
  available <- intersect(all_markers, rownames(obj))
  
  validation_list <- list()
  
  for (annot in annotations) {
    cells <- which(obj@meta.data$cluster_annotation == annot)
    n_cells <- length(cells)
    
    if (n_cells < 10) next
    
    # Find matching markers from config
    pos_markers <- c()
    neg_markers <- c()
    
    # Search all categories
    for (cat in c("lineages", "cd4_subtypes", "cd8_subtypes", "nk_subtypes")) {
      if (cat %in% names(markers_config)) {
        for (ct in names(markers_config[[cat]])) {
          if (ct == annot || harmonization[[ct]] == annot) {
            pos_markers <- c(pos_markers, markers_config[[cat]][[ct]]$positive)
            neg_markers <- c(neg_markers, markers_config[[cat]][[ct]]$negative)
          }
        }
      }
    }
    
    pos_markers <- intersect(unique(pos_markers), available)
    neg_markers <- intersect(unique(neg_markers), available)
    
    # Calculate validation scores
    pos_score <- NA
    neg_score <- NA
    
    if (length(pos_markers) > 0) {
      pos_expr <- sapply(pos_markers, function(m) {
        mean(obj@meta.data[[paste0("Score_", annot)]][cells] > 0.3, na.rm = TRUE)
      })
      pos_score <- mean(pos_expr, na.rm = TRUE)
    }
    
    validation_list[[annot]] <- data.frame(
      annotation = annot,
      n_cells = n_cells,
      n_pos_markers = length(pos_markers),
      n_neg_markers = length(neg_markers),
      validation_score = ifelse(!is.na(pos_score), pos_score, 0),
      stringsAsFactors = FALSE
    )
  }
  
  validation_df <- do.call(rbind, validation_list)
  rownames(validation_df) <- NULL
  
  log_message("Validation complete")
  validation_df
}