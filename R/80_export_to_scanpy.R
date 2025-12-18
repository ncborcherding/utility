# R/80_export_to_scanpy.R
# =============================================================================
# Export Seurat object to H5AD format for Python/scanpy users
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(yaml)
  library(BPCells)
  library(reticulate)
})

source("R/00_utils.R")

# =============================================================================
# SETUP CONDA ENVIRONMENT
# =============================================================================

setup_python_env <- function() {
  tryCatch({
    use_condaenv("sc-integration-benchmark", required = TRUE)
    
    # Verify required packages
    if (!py_module_available("h5py") || !py_module_available("numpy")) {
      stop("Python packages h5py and numpy are required")
    }
    
    log_message("Python environment ready")
    return(TRUE)
  }, error = function(e) {
    stop("Failed to set up Python environment: ", e$message)
  })
}


# =============================================================================
# PYTHON-BASED H5AD WRITER (handles int64 properly)
# =============================================================================

#' Initialize an empty H5AD file with proper structure
init_h5ad_file <- function(h5ad_path, n_cells, n_genes) {
  
  log_message("  Initializing H5AD structure...")
  
  py_run_string(sprintf('
import h5py
import numpy as np

h5ad_path = "%s"
n_cells = %d
n_genes = %d

with h5py.File(h5ad_path, "w") as f:
    # Create X group for CSR sparse matrix
    X = f.create_group("X")
    X.attrs["encoding-type"] = "csr_matrix"
    X.attrs["encoding-version"] = "0.1.0"
    X.attrs.create("shape", data=np.array([n_cells, n_genes], dtype=np.int64))
    
    # Create datasets with chunking and compression
    # Use int64 for indices and indptr to handle >2^31 elements
    X.create_dataset("data", shape=(0,), maxshape=(None,), 
                     dtype=np.float32, chunks=(1000000,), 
                     compression="gzip", compression_opts=4)
    X.create_dataset("indices", shape=(0,), maxshape=(None,), 
                     dtype=np.int64, chunks=(1000000,), 
                     compression="gzip", compression_opts=4)
    X.create_dataset("indptr", shape=(n_cells + 1,), 
                     dtype=np.int64, compression="gzip", compression_opts=4)
    
    # Initialize indptr[0] = 0
    f["X/indptr"][0] = 0
    
    # Create other required groups
    for grp in ["obs", "var", "obsm", "varm", "obsp", "varp", "uns", "layers"]:
        g = f.create_group(grp)
        if grp in ["obsm", "varm", "obsp", "varp", "uns", "layers"]:
            g.attrs["encoding-type"] = "dict"
            g.attrs["encoding-version"] = "0.1.0"

print("  H5AD structure initialized")
', h5ad_path, n_cells, n_genes))
}


#' Append a chunk of sparse data using Python (int64 safe)
append_sparse_chunk <- function(h5ad_path, data_vec, indices_vec, 
                                indptr_cumsum, start_row, end_row, 
                                chunk_num, total_chunks) {
  
  # Handle empty chunks
  if (length(data_vec) == 0) {
    log_message(sprintf("    Chunk %d/%d: empty (no non-zeros)", chunk_num, total_chunks))
    return(invisible(NULL))
  }
  
  # Ensure vectors are plain numeric (not Matrix class objects)
  data_vec <- as.numeric(data_vec)
  indices_vec <- as.numeric(indices_vec)
  indptr_cumsum <- as.numeric(indptr_cumsum)
  
  # Save arrays to temp binary files (avoids reticulate memory issues)
  temp_dir <- tempdir()
  data_file <- file.path(temp_dir, "chunk_data.bin")
  indices_file <- file.path(temp_dir, "chunk_indices.bin")
  indptr_file <- file.path(temp_dir, "chunk_indptr.bin")
  
  # Write as binary: float32 for data, float64 for indices (convert to int64 in Python)
  con <- file(data_file, "wb")
  writeBin(data_vec, con, size = 4)  # float32
  close(con)
  
  con <- file(indices_file, "wb")
  writeBin(indices_vec, con, size = 8)  # float64 -> int64 in Python
  close(con)
  
  con <- file(indptr_file, "wb")
  writeBin(indptr_cumsum, con, size = 8)  # float64 -> int64 in Python
  close(con)
  
  py_run_string(sprintf('
import h5py
import numpy as np

# Read from temp files
data = np.fromfile("%s", dtype=np.float32)
indices = np.fromfile("%s", dtype=np.float64).astype(np.int64)
indptr_cumsum = np.fromfile("%s", dtype=np.float64).astype(np.int64)

h5ad_path = "%s"
start_row = %d
end_row = %d
chunk_num = %d
total_chunks = %d

with h5py.File(h5ad_path, "a") as f:
    # Current sizes
    curr_nnz = f["X/data"].shape[0]
    
    # Resize and append data
    new_nnz = curr_nnz + len(data)
    f["X/data"].resize((new_nnz,))
    f["X/data"][curr_nnz:new_nnz] = data
    
    # Resize and append indices
    f["X/indices"].resize((new_nnz,))
    f["X/indices"][curr_nnz:new_nnz] = indices
    
    # Write indptr for this chunk (offset by current nnz)
    # indptr[i] = cumulative nnz up through row i-1
    # For rows start_row..end_row, we fill indptr[start_row+1..end_row+1]
    indptr_offset = indptr_cumsum + curr_nnz
    f["X/indptr"][start_row + 1 : end_row + 2] = indptr_offset

if chunk_num %% 10 == 0 or chunk_num == total_chunks:
    print(f"    Chunk {chunk_num}/{total_chunks}: {new_nnz:,} total non-zeros")
',
                        gsub("\\\\", "/", data_file),
                        gsub("\\\\", "/", indices_file),
                        gsub("\\\\", "/", indptr_file),
                        h5ad_path, start_row, end_row, chunk_num, total_chunks
  ))
  
  # Cleanup
  unlink(c(data_file, indices_file, indptr_file))
}


#' Write obs (cell metadata) via Python using anndata-native approach
#' This delegates to anndata to ensure categorical encoding is compatible
write_obs_python <- function(h5ad_path, meta_data, cell_names) {
  
  log_message("  Writing cell metadata (via anndata)...")
  
  # Save to temp CSV
  temp_csv <- tempfile(fileext = ".csv")
  meta_out <- meta_data
  meta_out[["_row_names_"]] <- cell_names
  fwrite(meta_out, temp_csv)
  
  py_run_string(sprintf('
import h5py
import numpy as np
import pandas as pd
import anndata as ad
import tempfile
import os
import warnings
warnings.filterwarnings("ignore")

# Read metadata
df = pd.read_csv("%s", low_memory=False)
df = df.set_index("_row_names_")
df.index.name = None

# Convert object columns to categorical (anndata will handle encoding)
for col in df.columns:
    if df[col].dtype == object:
        df[col] = df[col].fillna("NA").astype("category")

# Create a minimal AnnData just for obs
n_cells = len(df)
dummy_X = np.zeros((n_cells, 1), dtype=np.float32)
adata_temp = ad.AnnData(X=dummy_X, obs=df)

h5ad_path = "%s"

# Remove the empty obs group we created during init
with h5py.File(h5ad_path, "r+") as f:
    if "obs" in f:
        del f["obs"]

# Write temp AnnData, then copy just the obs group
temp_h5ad = tempfile.mktemp(suffix=".h5ad")
adata_temp.write_h5ad(temp_h5ad)

# Copy obs from temp file to our main H5AD
with h5py.File(temp_h5ad, "r") as src:
    with h5py.File(h5ad_path, "r+") as dst:
        src.copy("obs", dst)

os.remove(temp_h5ad)

print(f"  Wrote {len(df.columns)} metadata columns")
', gsub("\\\\", "/", temp_csv), h5ad_path))
  
  unlink(temp_csv)
}


#' Write var (gene names) via Python
write_var_python <- function(h5ad_path, gene_names) {
  
  log_message("  Writing gene names...")
  
  temp_file <- tempfile(fileext = ".txt")
  writeLines(gene_names, temp_file)
  
  py_run_string(sprintf('
import h5py
import numpy as np

with open("%s", "r") as f:
    gene_names = [line.strip() for line in f]

with h5py.File("%s", "a") as f:
    var = f["var"]
    var.create_dataset("_index", data=np.array(gene_names, dtype="S"), compression="gzip")
    var.attrs["encoding-type"] = "dataframe"
    var.attrs["encoding-version"] = "0.2.0"
    var.attrs["_index"] = "_index"
    var.attrs.create("column-order", data=np.array([], dtype="S"))

print(f"  Wrote {len(gene_names):,} gene names")
', gsub("\\\\", "/", temp_file), h5ad_path))
  
  unlink(temp_file)
}


#' Write obsm (embeddings) via Python
write_obsm_python <- function(h5ad_path, obj) {
  
  log_message("  Writing embeddings...")
  
  n_cells <- ncol(obj)
  
  for (red_name in names(obj@reductions)) {
    emb <- Embeddings(obj, reduction = red_name)
    if (!is.matrix(emb)) emb <- as.matrix(emb)
    if (nrow(emb) != n_cells) next
    
    # Standard naming for anndata
    h5_name <- switch(tolower(red_name),
                      "pca" = "X_pca",
                      "umap" = "X_umap", 
                      "harmony" = "X_harmony",
                      "scvi" = "X_scvi",
                      "scanvi" = "X_scanvi",
                      "fastmnn" = "X_fastmnn",
                      paste0("X_", red_name))
    
    # Write via temp file - convert to row-major numeric vector
    temp_file <- tempfile(fileext = ".bin")
    emb_vec <- as.numeric(as.vector(t(emb)))  # Row-major order for numpy
    
    con <- file(temp_file, "wb")
    writeBin(emb_vec, con, size = 4)  # float32
    close(con)
    
    tryCatch({
      py_run_string(sprintf('
import h5py
import numpy as np

emb = np.fromfile("%s", dtype=np.float32).reshape(%d, %d)

with h5py.File("%s", "a") as f:
    f["obsm"].create_dataset("%s", data=emb, compression="gzip", compression_opts=4)

print(f"    %s: {emb.shape}")
', gsub("\\\\", "/", temp_file), nrow(emb), ncol(emb), h5ad_path, h5_name, h5_name))
    }, error = function(e) {
      warning(sprintf("Could not write %s: %s", red_name, e$message))
    })
    
    unlink(temp_file)
  }
}


# =============================================================================
# VERIFICATION
# =============================================================================

#' Verify H5AD file loads correctly with anndata
verify_h5ad <- function(h5ad_path) {
  
  tryCatch({
    py_run_string(sprintf('
import anndata as ad
import warnings
warnings.filterwarnings("ignore")

h5ad_path = r"%s"
verification_passed = False

try:
    adata = ad.read_h5ad(h5ad_path)
    print(f"  ✓ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Check sparse matrix
    if hasattr(adata.X, "indptr"):
        if adata.X.indptr.min() < 0:
            print("  ✗ Integer overflow in indptr!")
        else:
            print(f"  ✓ Sparse matrix: {adata.X.nnz:,} non-zeros")
    
    # Check embeddings
    print(f"  ✓ Embeddings: {list(adata.obsm.keys())}") 
    
    verification_passed = True
    
except Exception as e:
    print(f"  ✗ Failed: {e}")
    verification_passed = False
', h5ad_path))
    
    # Get result from Python
    result <- py$verification_passed
    return(isTRUE(result))
    
  }, error = function(e) {
    log_message(sprintf("  ✗ Verification error: %s", e$message))
    return(FALSE)
  })
}


# =============================================================================
# MAIN EXPORT FUNCTION
# =============================================================================

#' Export Seurat object to H5AD
#' 
#' @param config_path Path to config.yaml
#' @param input_path Path to Seurat object (optional)
#' @param chunk_size Cells per chunk (default: 50000)
#' @param output_dir Output directory (optional)
#' @param create_h5mu Whether to create H5MU file (default: TRUE)
export_to_scanpy <- function(config_path = "config.yaml",
                             input_path = NULL,
                             chunk_size = 50000,
                             output_dir = NULL,
                             create_h5mu = TRUE) {
  
  log_message("═══════════════════════════════════════════════════════════════")
  log_message("  EXPORT TO H5AD (int64 indices, gzip compression)")
  log_message("═══════════════════════════════════════════════════════════════")
  
  # Setup Python
  setup_python_env()
  
  # Load config
  config <- yaml::read_yaml(config_path)
  
  if (is.null(input_path)) {
    input_path <- file.path(config$paths$results_dir, "FINAL_integrated_object.rds")
  }
  
  if (is.null(output_dir)) {
    output_dir <- file.path(config$paths$results_dir, "09_scanpy_export")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  h5ad_path <- file.path(output_dir, "adata.h5ad")
  
  # Load object
  log_message("Loading: ", basename(input_path))
  obj <- readRDS(input_path)
  obj$Frequency <- NULL
  obj$clonoType <- NULL
  
  if (isTRUE(obj@misc$lite_version)) {
    stop("Cannot export lite object. Use full BPCells version.")
  }
  
  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  cell_names <- colnames(obj)
  gene_names <- rownames(obj)
  
  log_message(sprintf("  %s cells × %s genes", 
                      format(n_cells, big.mark = ","),
                      format(n_genes, big.mark = ",")))
  
  # Get counts
  counts <- tryCatch(obj@assays$RNA@layers$counts, error = function(e) NULL)
  if (is.null(counts)) counts <- tryCatch(obj@assays$RNA@counts, error = function(e) NULL)
  if (is.null(counts)) stop("No count matrix found")
  
  is_bpcells <- inherits(counts, "IterableMatrix") || grepl("BPCells", class(counts)[1])
  log_message(sprintf("  Matrix: %s", if(is_bpcells) "BPCells" else class(counts)[1]))
  
  # Remove old file
  if (file.exists(h5ad_path)) unlink(h5ad_path)
  
  # Step 1: Initialize
  log_message("")
  log_message("[1/6] Creating H5AD structure...")
  init_h5ad_file(h5ad_path, n_cells, n_genes)
  
  # Step 2: Write sparse matrix in chunks
  log_message("")
  log_message("[2/6] Writing sparse matrix (chunked, int64)...")
  
  n_chunks <- ceiling(n_cells / chunk_size)
  log_message(sprintf("  %d chunks of ~%d cells", n_chunks, chunk_size))
  
  total_nnz <- 0
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_cells)
    
    # Extract chunk
    chunk_counts <- counts[, start_idx:end_idx, drop = FALSE]
    
    # Convert to dgCMatrix
    if (is_bpcells || !inherits(chunk_counts, "dgCMatrix")) {
      chunk_sparse <- as(chunk_counts, "dgCMatrix")
    } else {
      chunk_sparse <- chunk_counts
    }
    
    # Extract components directly - no transpose needed!
    data_vec <- as.numeric(chunk_sparse@x)
    indices_vec <- as.numeric(chunk_sparse@i)  # Gene indices (0-based)
    indptr_vec <- as.numeric(chunk_sparse@p)   # Cell boundaries
    
    # Calculate per-cell nnz for this chunk
    indptr_diff <- diff(indptr_vec)  # nnz per cell
    indptr_cumsum <- cumsum(indptr_diff)  # Cumulative for this chunk
    
    # Handle empty chunks
    if (length(data_vec) == 0) {
      data_vec <- numeric(0)
      indices_vec <- numeric(0)
    }
    
    chunk_nnz <- length(data_vec)
    total_nnz <- total_nnz + chunk_nnz
    
    # Write via Python
    append_sparse_chunk(h5ad_path, data_vec, indices_vec, indptr_cumsum,
                        start_idx - 1, end_idx - 1, i, n_chunks)
    
    rm(chunk_sparse, chunk_counts)
    if (i %% 5 == 0) gc()
  }
  
  log_message(sprintf("  Total: %s non-zeros", format(total_nnz, big.mark = ",")))
  
  # Step 3: Metadata
  log_message("")
  log_message("[3/6] Writing metadata & embeddings...")
  write_obs_python(h5ad_path, obj@meta.data, cell_names)
  write_var_python(h5ad_path, gene_names)
  write_obsm_python(h5ad_path, obj)
  
  # Step 4: Verify H5AD
  log_message("")
  log_message("[4/6] Verifying H5AD...")
  h5ad_ok <- verify_h5ad(h5ad_path)
  
  if (!h5ad_ok) {
    log_message("  ⚠ H5AD verification failed - will attempt repair during H5MU creation")
  }
  
  # Step 5: AIRR
  log_message("")
  log_message("[5/6] Exporting TCR data...")
  airr_path <- export_airr_data(obj, output_dir)
  
  # Step 6: Create H5MU (optional)
  h5mu_path <- NULL
  if (create_h5mu && !is.null(airr_path)) {
    log_message("")
    log_message("[6/6] Creating H5MU (MuData)...")
    h5mu_path <- create_h5mu_with_reticulate(h5ad_path, airr_path, output_dir)
    
    if (is.null(h5mu_path)) {
      log_message("  H5MU creation failed")
    }
  } else if (create_h5mu && is.null(airr_path)) {
    log_message("")
    log_message("[6/6] Skipping H5MU (no AIRR data)")
  } else {
    log_message("")
    log_message("[6/6] Skipping H5MU (create_h5mu=FALSE)")
  }
  
  # Generate helpers
  generate_verification_script(output_dir)
  generate_info(output_dir, n_cells, n_genes, !is.null(airr_path), total_nnz)
  
  # Final verification
  log_message("")
  log_message("Running final verification...")
  final_ok <- verify_h5ad(h5ad_path)
  
  # Summary
  file_size_gb <- file.size(h5ad_path) / 1e9
  h5mu_size_gb <- if (!is.null(h5mu_path) && file.exists(h5mu_path)) file.size(h5mu_path) / 1e9 else NA
  
  log_message("")
  log_message("═══════════════════════════════════════════════════════════════")
  log_message("  EXPORT COMPLETE")
  log_message("═══════════════════════════════════════════════════════════════")
  log_message("")
  log_message("  Files created:")
  log_message(sprintf("    adata.h5ad               %.2f GB %s", file_size_gb, if(final_ok) "✓" else "⚠"))
  if (!is.null(airr_path)) {
    airr_size_mb <- file.size(airr_path) / 1e6
    log_message(sprintf("    airr_rearrangement.tsv.gz  %.1f MB", airr_size_mb))
  }
  if (!is.null(h5mu_path) && file.exists(h5mu_path)) {
    log_message(sprintf("    mudata.h5mu              %.2f GB", h5mu_size_gb))
  }
  log_message("")
  log_message(sprintf("  Non-zeros: %s (int64 indices)", format(total_nnz, big.mark = ",")))
  
  return(invisible(list(
    h5ad = h5ad_path,
    airr = airr_path,
    h5mu = h5mu_path
  )))
}


# =============================================================================
# AIRR EXPORT
# =============================================================================

export_airr_data <- function(obj, output_dir) {
  
  meta <- obj@meta.data
  
  if (!"CTgene" %in% colnames(meta)) {
    log_message("  No TCR data found")
    return(NULL)
  }
  
  log_message("  Parsing TCR data...")
  
  airr_list <- vector("list", nrow(meta))
  valid_count <- 0
  
  for (i in seq_len(nrow(meta))) {
    ct_gene <- meta$CTgene[i]
    if (is.na(ct_gene) || ct_gene == "" || ct_gene == "NA") next
    
    ct_aa <- if("CTaa" %in% colnames(meta)) meta$CTaa[i] else NA
    ct_nt <- if("CTnt" %in% colnames(meta)) meta$CTnt[i] else NA
    
    gene_parts <- strsplit(as.character(ct_gene), "_")[[1]]
    aa_parts <- if(!is.na(ct_aa) && ct_aa != "") strsplit(as.character(ct_aa), "_")[[1]] else character(0)
    nt_parts <- if(!is.na(ct_nt) && ct_nt != "") strsplit(as.character(ct_nt), "_")[[1]] else character(0)
    
    cell_id <- rownames(meta)[i]
    
    for (j in seq_along(gene_parts)) {
      gene_str <- gene_parts[j]
      if (gene_str == "" || gene_str == "NA") next
      
      genes <- strsplit(gene_str, "\\.")[[1]]
      v_gene <- genes[1]
      
      locus <- if(grepl("^TRA", v_gene)) "TRA" else
        if(grepl("^TRB", v_gene)) "TRB" else
          if(grepl("^TRG", v_gene)) "TRG" else
            if(grepl("^TRD", v_gene)) "TRD" else
              if(grepl("^IGH", v_gene)) "IGH" else
                if(grepl("^IGK", v_gene)) "IGK" else
                  if(grepl("^IGL", v_gene)) "IGL" else NA
      
      if (is.na(locus)) next
      
      valid_count <- valid_count + 1
      airr_list[[valid_count]] <- list(
        cell_id = cell_id,
        sequence_id = paste0(cell_id, "_", j),
        locus = locus,
        productive = TRUE,
        v_call = v_gene,
        j_call = if(length(genes) > 1) genes[length(genes)] else NA_character_,
        junction_aa = if(j <= length(aa_parts)) aa_parts[j] else NA_character_,
        junction = if(j <= length(nt_parts)) nt_parts[j] else NA_character_
      )
    }
  }
  
  if (valid_count == 0) {
    log_message("  No valid TCR chains")
    return(NULL)
  }
  
  airr_df <- rbindlist(airr_list[1:valid_count], fill = TRUE)
  
  # Save as gzip-compressed TSV (fwrite handles .gz extension automatically)
  airr_path <- file.path(output_dir, "airr_rearrangement.tsv.gz")
  fwrite(airr_df, airr_path, sep = "\t", compress = "gzip")
  
  # Report file size
  
  file_size_mb <- file.size(airr_path) / 1e6
  log_message(sprintf("  %s chains from %s cells (%.1f MB compressed)",
                      format(nrow(airr_df), big.mark = ","),
                      format(length(unique(airr_df$cell_id)), big.mark = ","),
                      file_size_mb))
  
  return(airr_path)
}


# =============================================================================
# HELPER GENERATORS
# =============================================================================

generate_verification_script <- function(output_dir) {
  script <- '#!/usr/bin/env python3
"""Verify H5AD export with int64 index check."""
import sys
try:
    import scanpy as sc
    import numpy as np
except ImportError as e:
    print(f"Missing: {e}")
    sys.exit(1)

print(f"scanpy: {sc.__version__}")
print("\\nLoading adata.h5ad...")

try:
    adata = sc.read_h5ad("adata.h5ad")
    print(f"✓ Shape: {adata.n_obs:,} × {adata.n_vars:,}")
    print(f"  obsm: {list(adata.obsm.keys())}")
    
    if hasattr(adata.X, "indptr"):
        print(f"  indptr dtype: {adata.X.indptr.dtype}")
        print(f"  indices dtype: {adata.X.indices.dtype}")
        print(f"  nnz: {adata.X.nnz:,}")
        
        if adata.X.indptr.min() < 0:
            print("✗ ERROR: Integer overflow detected!")
            sys.exit(1)
        print("✓ No overflow - indices are valid")
    
    print("\\n✓ Verification passed!")
except Exception as e:
    print(f"✗ {e}")
    sys.exit(1)
'
  writeLines(script, file.path(output_dir, "verify_export.py"))
}

generate_readme <- function(output_dir, n_cells, n_genes, has_airr, total_nnz) {
  readme <- sprintf('uTILity Export
==============
%s cells × %s genes
%s non-zeros (int64 indices)

Files:
  adata.h5ad - Gene expression
  %s
  verify_export.py - Verification

Python:
  import scanpy as sc
  adata = sc.read_h5ad("adata.h5ad")

R (zellkonverter):
  library(zellkonverter)
  sce <- readH5AD("adata.h5ad")
  obj <- as.Seurat(sce, counts="X", data=NULL)

R (SeuratDisk):
  library(SeuratDisk)
  Convert("adata.h5ad", dest="h5seurat")
  obj <- LoadH5Seurat("adata.h5seurat")
%s
Verify: python verify_export.py
',
                    format(n_cells, big.mark = ","),
                    format(n_genes, big.mark = ","),
                    format(total_nnz, big.mark = ","),
                    if(has_airr) "airr_rearrangement.tsv.gz - TCR (AIRR format, gzip)" else "",
                    if(has_airr) "\nTCR (Python):\n  import pandas as pd\n  import scirpy as ir\n  # Load gzipped AIRR via pandas\n  airr_df = pd.read_csv('airr_rearrangement.tsv.gz', sep='\\\\t', compression='gzip')\n  airr_df.to_csv('temp.tsv', sep='\\\\t', index=False)\n  adata_tcr = ir.io.read_airr('temp.tsv')" else ""
  )
  writeLines(readme, file.path(output_dir, "INFO.txt"))
  
  if (has_airr) {
    h5mu_script <- '#!/usr/bin/env python3
"""Create H5MU (optional - duplicates data, large file)"""
import scanpy as sc
import scirpy as ir
import muon as mu
import pandas as pd
import tempfile
import os

print("Loading GEX...")
adata_gex = sc.read_h5ad("adata.h5ad")

print("Loading AIRR...")
# scirpy may have issues with gzipped files, so read via pandas first
airr_file = "airr_rearrangement.tsv.gz"
if airr_file.endswith(".gz"):
    airr_df = pd.read_csv(airr_file, sep="\\t", compression="gzip")
    temp_tsv = tempfile.mktemp(suffix=".tsv")
    airr_df.to_csv(temp_tsv, sep="\\t", index=False)
    adata_airr = ir.io.read_airr(temp_tsv)
    os.remove(temp_tsv)
else:
    adata_airr = ir.io.read_airr(airr_file)

print("Creating MuData...")
mdata = mu.MuData({"gex": adata_gex, "airr": adata_airr})
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

print("Saving...")
mdata.write("mudata.h5mu")
print("Done: mudata.h5mu")
'
    writeLines(h5mu_script, file.path(output_dir, "create_h5mu.py"))
  }
}


#' Create H5MU combining GEX and AIRR data
create_h5mu_with_reticulate <- function(h5ad_path, airr_path, export_dir) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    log_message("  reticulate not available")
    return(NULL)
  }
  
  # Check required packages
  required_pkgs <- c("anndata", "muon", "scirpy", "h5py", "numpy")
  for (pkg in required_pkgs) {
    if (!reticulate::py_module_available(pkg)) {
      log_message(sprintf("  Missing Python package: %s", pkg))
      return(NULL)
    }
  }
  
  h5mu_path <- file.path(export_dir, "mudata.h5mu")
  
  tryCatch({
    log_message("  Loading and creating MuData...")
    
    py_run_string(sprintf('
import anndata as ad
import scirpy as ir
from muon import MuData
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

print("  Loading GEX...")
adata_gex = ad.read_h5ad(r"%s")
print(f"  GEX: {adata_gex.n_obs:,} cells x {adata_gex.n_vars:,} genes")

print("  Loading AIRR...")
# Read gzipped TSV via pandas first, then pass to scirpy
airr_path = r"%s"
if airr_path.endswith(".gz"):
    airr_df = pd.read_csv(airr_path, sep="\\t", compression="gzip")
    # Save to temp uncompressed file for scirpy
    import tempfile
    temp_tsv = tempfile.mktemp(suffix=".tsv")
    airr_df.to_csv(temp_tsv, sep="\\t", index=False)
    adata_airr = ir.io.read_airr(temp_tsv)
    import os
    os.remove(temp_tsv)
else:
    adata_airr = ir.io.read_airr(airr_path)
print(f"  AIRR: {adata_airr.n_obs:,} cells")

print("  Creating MuData...")
mdata = MuData({"gex": adata_gex, "airr": adata_airr})

print("  Running scirpy...")
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

print("  Saving H5MU (gzip)...")
mdata.write(r"%s", compression="gzip")
print(f"  ✓ Saved: %s")
', h5ad_path, airr_path, h5mu_path, basename(h5mu_path)))
    
    log_message("  H5MU created successfully")
    return(h5mu_path)
    
  }, error = function(e) {
    log_message(sprintf("  Error creating H5MU: %s", e$message))
    return(NULL)
  })
}

# =============================================================================
if (sys.nframe() == 0) {
  log_message("Usage: export_to_scanpy('config.yaml')")
}