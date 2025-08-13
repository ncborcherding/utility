# install_dependencies.R
#
# This script installs all the necessary R packages for the uTILity analysis pipeline.
# It handles packages from CRAN and GitHub.

# --- Setup ---

# Create a user-specific library path if it doesn't exist.
# This avoids permission issues in shared environments.
lib_path <- "~/R_libs"
if (!dir.exists(lib_path)) {
  dir.create(lib_path, recursive = TRUE)
  cat("Created personal library at:", lib_path, "\n")
}
.libPaths(c(lib_path, .libPaths()))
cat("Using library paths:", paste(.libPaths(), collapse = "\n"), "\n")

# Increase timeout for slow-compiling packages
options(timeout = 3600)

# --- Package Installation ---

# Helper function to install packages if they are not already present
install_if_missing <- function(pkg_name, repo_type, repo_path = NULL) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat("Installing", pkg_name, "...\n")
    tryCatch({
      switch(repo_type,
        "cran" = install.packages(pkg_name, repos = "http://cran.us.r-project.org"),
        "github" = remotes::install_github(repo_path)
      )
      cat(pkg_name, "installed successfully.\n")
    }, error = function(e) {
      cat("ERROR: Failed to install", pkg_name, ".\n")
      print(e)
    })
  } else {
    cat(pkg_name, "is already installed.\n")
  }
}

# Install 'remotes' first, as it's needed for GitHub installations
install_if_missing("remotes", "cran")
library(remotes)

# --- CRAN Packages ---
cat("\n--- Installing CRAN Packages ---\n")
cran_packages <- c(
  "Seurat",      # Core analysis toolkit
  "cluster",     # For calculating silhouette scores
  "reticulate"   # For Python integration (scvi-tools)
)
for (pkg in cran_packages) {
  install_if_missing(pkg, "cran")
}

# --- GitHub Packages ---
cat("\n--- Installing GitHub Packages ---\n")
# BPCells for memory-efficient on-disk matrices
install_if_missing("BPCells", "github", "bnprks/BPCells")

# SeuratWrappers for scVI integration
install_if_missing("SeuratWrappers", "github", "satijalab/seurat-wrappers")


cat("\nDependency installation script finished.\n")
