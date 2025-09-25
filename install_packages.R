# install_packages.R

# This script checks for and installs all required R packages for the pipeline.

# --- CRAN Packages ---
cran_packages <- c(
  "yaml",          # For reading config file
  "dplyr",         # General data manipulation
  "purrr",         # Functional programming
  "ggplot2",       # Plotting
  "patchwork",     # Combining ggplot plots
  "RColorBrewer",  # Color palettes
  "viridis",       # More color palettes
  "Matrix",        # Sparse matrices
  "FNN",           # Fast nearest neighbor search
  "mclust",        # For ARI metric
  "aricode",       # For NMI and other metrics
  "reticulate",    # R interface to Python
  "future",        # Parallel processing
  "future.apply",  # Parallel apply functions
  "pak",           # Alternative package manager
  "remotes",       # For installing from GitHub
  "rjson"
)

# --- Bioconductor Packages ---
bioc_packages <- c(
  "Seurat",        # Core single-cell analysis
  "harmony",       # Harmony integration
  "batchelor",     # fastMNN integration
  "SingleCellExperiment", # Data structure for batchelor/scran
  "scran",         # HVG detection
  "scDblFinder",   # Doublet detection 
  "SingleR",       # Cell type annotation 
  "celldex",       # Annotation references 
  "Azimuth",       # Annotation 
  "SeuratData",    # For Azimuth references
  "BiocParallel",  # Parallel processing for Bioconductor
  "HGNChelper",     # For gene symbol checking
  "scRepertoire"
)

# --- GitHub Packages ---
github_packages <- c(
  "mengxu98/LISI"               # For LISI metric
)

# --- Installation Logic ---

# Function to install from CRAN if not already installed
install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "from CRAN..."))
    install.packages(pkg, ask = FALSE)
  } else {
    message(paste(pkg, "is already installed."))
  }
}

# Function to install from Bioconductor if not already installed
install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "from Bioconductor..."))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE)
  } else {
    message(paste(pkg, "is already installed."))
  }
}

# Function to install from GitHub if not already installed
install_if_missing_github <- function(pkg) {
  pkg_name <- basename(pkg)
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message(paste("Installing", pkg, "from GitHub..."))
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github(pkg)
  } else {
    message(paste(pkg_name, "is already installed."))
  }
}

# --- Run Installations ---

message("--- Checking and installing CRAN packages ---")
sapply(cran_packages, install_if_missing_cran)

message("\n--- Checking and installing Bioconductor packages ---")
sapply(bioc_packages, install_if_missing_bioc)

message("\n--- Checking and installing GitHub packages ---")
sapply(github_packages, install_if_missing_github)


# --- Special Case for Seurat v5 Objects ---
# The original code uses `options(Seurat.object.assay.version = "v5")`
# This requires the SeuratObject package. Seurat > 5.0 should install it.
# We also need SeuratDisk for converting to h5ad.
message("\n--- Checking for Seurat v5 dependencies ---")
install_if_missing_cran("SeuratDisk")


message("\n\nAll required packages have been checked and installed.")
