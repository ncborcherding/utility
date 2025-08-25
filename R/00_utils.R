# R/00_utils.R
# Helper functions for the integration benchmark pipeline

# --- Logging and Timing ---

#' Log a message with a timestamp
#'
#' @param ... Items to be pasted together into the log message.
log_message <- function(...) {
  message(paste0("[", Sys.time(), "] ", ...))
}

#' Start a timer
#'
#' @return The start time.
start_timer <- function() {
  tic <- Sys.time()
  log_message("Timer started.")
  return(tic)
}

#' Stop a timer and print the elapsed time
#'
#' @param tic The start time object from start_timer().
#' @param task_name A descriptive name for the timed task.
stop_timer <- function(tic, task_name = "Task") {
  toc <- Sys.time()
  elapsed <- as.numeric(difftime(toc, tic, units = "secs"))
  log_message(paste0(task_name, " finished. Time elapsed: ", round(elapsed, 2), " seconds."))
}


# --- File I/O ---

#' Save an R object to an RDS file, creating the directory if it doesn't exist.
#'
#' @param object The R object to save.
#' @param file_path The full path to the output file.
safe_save_rds <- function(object, file_path) {
  dir_name <- dirname(file_path)
  if (!dir.exists(dir_name)) {
    log_message("Creating directory: ", dir_name)
    dir.create(dir_name, recursive = TRUE)
  }
  log_message("Saving object to: ", file_path)
  saveRDS(object, file_path)
}


# --- Factor Sanitization ---

#' Sanitize factor levels to be valid variable names
#'
#' Replaces special characters with dots.
#' @param f A factor or character vector.
#' @return A character vector with sanitized names.
sanitize_factors <- function(f) {
  gsub("[^a-zA-Z0-9_.-]", ".", as.character(f))
}


# --- Plotting Theme ---

#' A custom ggplot theme for the project
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @param grid_lines Character indicating which grid lines to show ('Y', 'X', 'XY', 'No').
#' @param axis_lines Boolean, whether to draw axis lines.
#' @param legend_position Position of the legend.
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
