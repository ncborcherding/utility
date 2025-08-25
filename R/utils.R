.sample_cells <- function(n, pool) {
  if (n >= length(pool)) return(pool)
  sample(pool, n, replace = FALSE)
}

.silhouette_mean <- function(labels, emb, max_n = 4000L) {
  # labels: named factor/character, names = cells present in emb
  cells <- intersect(names(labels), rownames(emb))
  if (length(cells) < 20L || length(unique(labels[cells])) < 2L) return(NA_real_)
  if (length(cells) > max_n) cells <- sample(cells, max_n)
  d <- stats::dist(emb[cells, , drop = FALSE])  # Euclidean on embedding
  lab_int <- as.integer(factor(labels[cells]))
  sil <- cluster::silhouette(lab_int, d)
  mean(sil[, "sil_width"])
}

run_leiden_on_subset <- function(obj, cells, resolution) {
  sobj <- subset(obj, cells = cells)
  # Build neighbors in subset only
  sobj <- FindNeighbors(
    sobj,
    reduction  = reduction_for_dist,
    dims       = dims_use,
    k.param    = k_param,
    nn.method  = nn_method,
    verbose    = FALSE
  )
  sobj <- FindClusters(
    sobj,
    resolution = resolution,
    algorithm  = 4,          # Leiden
    random.seed = random_seed,
    verbose    = FALSE
  )
  list(
    ids    = Idents(sobj),   # factor, names = subset cells
    object = sobj
  )
}

UtilityTheme <- function(base_size = 12,
                             base_family = "sans",
                             grid_lines = "Y",
                             axis_lines = FALSE,
                             legend_position = "right") {
  
  t <- ggplot2::theme_bw(base_size = base_size, base_family = base_family)
  t <- t %+replace%
    ggplot2::theme(
      # Plot titles and caption
      plot.title = ggplot2::element_text(
        size = rel(1.2), hjust = 0, face = "bold",
        margin = ggplot2::margin(b = base_size / 2)
      ),
      plot.subtitle = ggplot2::element_text(
        size = rel(1.0), hjust = 0,
        margin = ggplot2::margin(b = base_size)
      ),
      plot.caption = ggplot2::element_text(
        size = rel(0.8), hjust = 1, color = "grey50",
        margin = ggplot2::margin(t = base_size / 2)
      ),
      
      # Backgrounds and borders
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.75),
      
      # Remove all grid lines by default; they will be added back conditionally
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Axis text, titles, and ticks
      axis.title = ggplot2::element_text(size = rel(1.0)),
      axis.text = ggplot2::element_text(size = rel(0.9), color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5),
      
      # Legend customization
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = rel(0.9), face = "bold"),
      legend.text = ggplot2::element_text(size = rel(0.85)),
      legend.position = legend_position,
      
      # Facet (strip) customization
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black", linewidth = 0.75),
      strip.text = ggplot2::element_text(
        size = rel(1.0), face = "bold", color = "black",
        margin = ggplot2::margin(t = base_size / 4, b = base_size / 4)
      )
    )
  
  # Conditionally add major grid lines based on the 'grid_lines' parameter
  grid_lines <- toupper(grid_lines)
  if (grid_lines %in% c("Y", "XY")) {
    t <- t + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.5))
  }
  if (grid_lines %in% c("X", "XY")) {
    t <- t + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "grey85", linewidth = 0.5))
  }
  
  # Conditionally add axis lines
  if (axis_lines) {
    t <- t + ggplot2::theme(axis.line = ggplot2::element_line(color = "black", linewidth = 0.5))
  }
  
  return(t)
}
