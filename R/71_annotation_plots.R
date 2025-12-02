# 51_annotation_plots.R
# ==============================================================================
# Visualization Functions for T Cell Annotation
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(viridis)
})

# ==============================================================================
# COLOR PALETTES
# ==============================================================================

#' Get standardized color palette for T cell annotations
#'
#' @param annotation_type Type of annotations (lineage, subtype, full)
#' @return Named vector of colors
get_tcell_palette <- function(annotation_type = "full") {
  
  lineage_colors <- c(
    "CD4_T" = "#E41A1C",
    "CD8_T" = "#377EB8", 
    "NK" = "#4DAF4A",
    "NKT" = "#984EA3",
    "gdT" = "#FF7F00",
    "MAIT" = "#FFFF33",
    "dnT" = "#A65628",
    "Unknown" = "#999999"
  )
  
  subtype_colors <- c(
    # CD4 subtypes
    "Th1" = "#E41A1C",
    "Th2" = "#FF6B6B",
    "Th17" = "#C91D42",
    "Th1/Th17" = "#D35400",
    "Th22" = "#E74C3C",
    "Tfh" = "#9B59B6",
    "Treg" = "#F39C12",
    "Treg_naive" = "#F5B041",
    "Treg_memory" = "#D68910",
    "CD4_naive" = "#FADBD8",
    "CD4_Tcm" = "#F1948A",
    "CD4_Tem" = "#E74C3C",
    "CD4_CTL" = "#922B21",
    "CD4_unassigned" = "#FFCCCB",
    
    # CD8 subtypes
    "CD8_naive" = "#AED6F1",
    "CD8_Tcm" = "#5DADE2",
    "CD8_Tem" = "#2E86C1",
    "CD8_Temra" = "#1B4F72",
    "CD8_exhausted" = "#154360",
    "CD8_CTL" = "#21618C",
    "CD8_proliferating" = "#3498DB",
    "CD8_unassigned" = "#D4E6F1",
    
    # NK subtypes
    "NK" = "#4DAF4A",
    "NK_CD56bright" = "#82E0AA",
    "NK_CD56dim" = "#239B56",
    "NK_adaptive" = "#145A32",
    "NK_unassigned" = "#D5F5E3",
    
    # Other
    "MAIT" = "#F4D03F",
    "gdT" = "#E67E22",
    "gdT_Vd2" = "#D35400",
    "gdT_nonVd2" = "#BA4A00",
    "dnT" = "#8B4513",
    "NKT" = "#9B59B6",
    
    # States
    "T_proliferating" = "#00CED1",
    "T_activated" = "#FF1493",
    "ISG_high" = "#9400D3",
    "Stressed" = "#696969",
    
    "Unknown" = "#999999"
  )
  
  if (annotation_type == "lineage") {
    return(lineage_colors)
  } else if (annotation_type == "subtype" || annotation_type == "full") {
    return(subtype_colors)
  } else {
    return(c(lineage_colors, subtype_colors))
  }
}

# ==============================================================================
# MARKER EXPRESSION PLOTS
# ==============================================================================

#' Create stacked violin plot for key markers
#'
#' @param obj Seurat object
#' @param features Vector of features to plot
#' @param group.by Grouping variable
#' @param cols Color palette
#' @return ggplot object
plot_stacked_violin <- function(obj, features, group.by = "cluster_annotation",
                                 cols = NULL) {
  
  # Get available features
  available <- intersect(features, rownames(obj))
  
  if (length(available) == 0) {
    stop("No features found in object")
  }
  
  # Use VlnPlot with stack
  VlnPlot(obj, features = available, group.by = group.by, 
          stack = TRUE, flip = TRUE, cols = cols) +
    theme(legend.position = "none")
}

#' Create feature plot grid for markers
#'
#' @param obj Seurat object
#' @param markers_by_type Named list of markers per cell type
#' @param reduction Reduction to use
#' @param ncol Number of columns
#' @return ggplot object
plot_marker_feature_grid <- function(obj, markers_by_type, 
                                      reduction = "umap", ncol = 4) {
  
  all_markers <- unlist(markers_by_type)
  available <- intersect(all_markers, rownames(obj))
  
  FeaturePlot(obj, features = available, reduction = reduction,
              ncol = ncol, order = TRUE) &
    scale_color_viridis_c() &
    theme(plot.title = element_text(size = 10))
}

#' Create comparative dot plot across annotation sources
#'
#' @param obj Seurat object
#' @param features Features to plot
#' @param group_cols Different grouping columns to compare
#' @return patchwork object
plot_comparative_dotplot <- function(obj, features,
                                      group_cols = c("cluster_annotation", 
                                                     "predicted.celltype.l2")) {
  
  available <- intersect(features, rownames(obj))
  plots <- list()
  
  for (gc in group_cols) {
    if (gc %in% colnames(obj@meta.data)) {
      p <- DotPlot(obj, features = available, group.by = gc) +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8)) +
        labs(title = gc)
      plots[[gc]] <- p
    }
  }
  
  wrap_plots(plots, nrow = 1)
}

# ==============================================================================
# CLUSTER COMPOSITION PLOTS
# ==============================================================================

#' Plot annotation composition per cluster
#'
#' @param obj Seurat object
#' @param cluster_col Cluster column
#' @param annotation_col Annotation column
#' @return ggplot object
plot_cluster_composition <- function(obj, cluster_col = "leiden.cluster",
                                      annotation_col = "predicted.celltype.l2") {
  
  comp_df <- obj@meta.data %>%
    count(.data[[cluster_col]], .data[[annotation_col]]) %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup()
  
  ggplot(comp_df, aes(x = factor(.data[[cluster_col]]), 
                      y = fraction, 
                      fill = .data[[annotation_col]])) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    labs(x = "Cluster", y = "Fraction", 
         title = paste("Composition by", annotation_col))
}

#' Plot annotation agreement between sources
#'
#' @param obj Seurat object
#' @param annot_cols Annotation columns to compare
#' @return ggplot object
plot_annotation_agreement <- function(obj, 
                                       annot_cols = c("cluster_annotation", 
                                                      "predicted.celltype.l2")) {
  
  if (length(annot_cols) != 2) {
    stop("Exactly 2 annotation columns required")
  }
  
  if (!all(annot_cols %in% colnames(obj@meta.data))) {
    stop("Annotation columns not found in metadata")
  }
  
  # Create confusion-style matrix
  conf_df <- obj@meta.data %>%
    count(.data[[annot_cols[1]]], .data[[annot_cols[2]]]) %>%
    group_by(.data[[annot_cols[1]]]) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup()
  
  ggplot(conf_df, aes(x = .data[[annot_cols[1]]], 
                       y = .data[[annot_cols[2]]], 
                       fill = fraction)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    labs(title = "Annotation Agreement",
         fill = "Fraction")
}

# ==============================================================================
# DIFFERENTIAL EXPRESSION HELPERS
# ==============================================================================

#' Find markers for each annotation
#'
#' @param obj Seurat object
#' @param annotation_col Annotation column
#' @param only.pos Only return positive markers
#' @param min.pct Minimum percent expression
#' @param logfc.threshold Log fold change threshold
#' @return Data frame with markers
find_annotation_markers <- function(obj, annotation_col = "cluster_annotation",
                                     only.pos = TRUE, min.pct = 0.25,
                                     logfc.threshold = 0.5) {
  
  Idents(obj) <- annotation_col
  
  markers <- FindAllMarkers(obj, 
                            only.pos = only.pos,
                            min.pct = min.pct,
                            logfc.threshold = logfc.threshold)
  
  markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 20) %>%
    ungroup()
}

#' Plot top markers heatmap
#'
#' @param obj Seurat object
#' @param markers Marker data frame from FindAllMarkers
#' @param n_top Number of top markers per group
#' @param annotation_col Annotation column
#' @return Heatmap
plot_top_markers_heatmap <- function(obj, markers, n_top = 5,
                                      annotation_col = "cluster_annotation") {
  
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = n_top) %>%
    pull(gene) %>%
    unique()
  
  DoHeatmap(obj, features = top_markers, group.by = annotation_col,
            size = 3) +
    theme(axis.text.y = element_text(size = 6))
}

# ==============================================================================
# QC VISUALIZATION
# ==============================================================================

#' Plot annotation confidence distribution
#'
#' @param obj Seurat object with annotation_confidence column
#' @param group.by Optional grouping
#' @return ggplot object
plot_confidence_distribution <- function(obj, group.by = NULL) {
  
  if (!"annotation_confidence" %in% colnames(obj@meta.data)) {
    stop("annotation_confidence column not found")
  }
  
  df <- obj@meta.data
  
  if (is.null(group.by)) {
    p <- ggplot(df, aes(x = annotation_confidence)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      labs(x = "Annotation Confidence", y = "Cell Count",
           title = "Distribution of Annotation Confidence")
  } else {
    p <- ggplot(df, aes(x = annotation_confidence, fill = .data[[group.by]])) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(x = "Annotation Confidence", y = "Density",
           title = "Annotation Confidence by Group")
  }
  
  p
}

#' Create summary figure with multiple panels
#'
#' @param obj Seurat object
#' @param markers_config Loaded marker configuration
#' @param output_file Output file path
#' @return Invisible NULL, saves plot
create_summary_figure <- function(obj, markers_config = NULL, 
                                   output_file = "annotation_summary.pdf") {
  
  # Panel A: UMAP with cluster annotation
  p1 <- DimPlot(obj, group.by = "cluster_annotation", 
                label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("Final Annotations") +
    theme(legend.position = "none")
  
  # Panel B: UMAP with leiden clusters
  p2 <- DimPlot(obj, group.by = "leiden.cluster",
                label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("Leiden Clusters") +
    theme(legend.position = "none")
  
  # Panel C: Key markers
  key_markers <- c("CD4", "CD8A", "FOXP3", "TBX21", "GATA3", "RORC",
                   "GZMB", "PRF1", "CCR7", "SELL")
  available <- intersect(key_markers, rownames(obj))
  
  p3 <- FeaturePlot(obj, features = available[1:min(4, length(available))],
                    ncol = 2, order = TRUE) &
    scale_color_viridis_c() &
    theme(plot.title = element_text(size = 8))
  
  # Panel D: Composition bar plot
  if ("predicted.celltype.l2" %in% colnames(obj@meta.data)) {
    p4 <- plot_cluster_composition(obj, "leiden.cluster", "predicted.celltype.l2") +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 6))
  } else {
    p4 <- ggplot() + theme_void() + labs(title = "No reference annotations")
  }
  
  # Combine
  combined <- (p1 | p2) / (p3 | p4)
  
  ggsave(output_file, combined, width = 16, height = 14)
  log_message("Summary figure saved to: ", output_file)
  
  invisible(NULL)
}

# ==============================================================================
# INTERACTIVE REPORTING
# ==============================================================================

#' Generate cluster annotation table for export
#'
#' @param obj Seurat object
#' @param cluster_col Cluster column
#' @param annotation_cols Annotation columns to include
#' @return Data frame summary
generate_cluster_table <- function(obj, cluster_col = "leiden.cluster",
                                    annotation_cols = c("cluster_annotation",
                                                        "predicted.celltype.l2",
                                                        "Monaco.labels")) {
  
  clusters <- sort(unique(obj@meta.data[[cluster_col]]))
  
  result <- data.frame(cluster = clusters)
  result$n_cells <- as.numeric(table(obj@meta.data[[cluster_col]])[as.character(clusters)])
  
  for (ac in annotation_cols) {
    if (ac %in% colnames(obj@meta.data)) {
      # Get majority for each cluster
      majority <- obj@meta.data %>%
        group_by(.data[[cluster_col]]) %>%
        count(.data[[ac]]) %>%
        slice_max(n, n = 1) %>%
        mutate(pct = round(n / sum(n) * 100, 1))
      
      result[[paste0(ac, "_majority")]] <- 
        majority[[ac]][match(result$cluster, majority[[cluster_col]])]
      result[[paste0(ac, "_pct")]] <- 
        majority$pct[match(result$cluster, majority[[cluster_col]])]
    }
  }
  
  result
}

#' Log message with timestamp
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
}
