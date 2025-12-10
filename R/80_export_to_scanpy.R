# R/80_export_to_scanpy.R
# 
# Unified export for both R and Python users
# Handles massive matrices (>2^31 non-zeros) via chunked export
#

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(yaml)
  library(BPCells)
  library(rhdf5)  # For direct HDF5 writing
})

source("R/00_utils.R")

# =============================================================================
# LOW-LEVEL HDF5 HELPERS (must be defined first)
# =============================================================================

#' Write a fixed-size integer array attribute (required for 'shape')
#' Uses h5py via reticulate to ensure proper INT64 type that anndata expects
write_int_array_attribute <- function(h5file, group_path, attr_name, values) {
  
  # First try using rhdf5's low-level API with explicit namespace
  tryCatch({
    fid <- H5Fopen(h5file)
    gid <- H5Gopen(fid, group_path)
    
    # Try to delete existing attribute if present
    if (H5Aexists(gid, attr_name)) {
      H5Adelete(gid, attr_name)
    }
    
    # Try low-level API for proper int64 type
    sid <- rhdf5::H5Screate_simple(dims = length(values))
    tid <- rhdf5::H5Tcopy("H5T_NATIVE_INT64")
    aid <- rhdf5::H5Acreate(gid, attr_name, tid, sid)
    rhdf5::H5Awrite(aid, as.integer(values))
    
    rhdf5::H5Aclose(aid)
    rhdf5::H5Tclose(tid)
    rhdf5::H5Sclose(sid)
    H5Gclose(gid)
    H5Fclose(fid)
    
    return(invisible(TRUE))
    
  }, error = function(e) {
    # Low-level API not available, try h5writeAttribute
    tryCatch({
      H5Gclose(gid)
      H5Fclose(fid)
    }, error = function(e2) {})
    
    # Fall back to h5writeAttribute - may need Python fix later
    tryCatch({
      fid <- H5Fopen(h5file)
      gid <- H5Gopen(fid, group_path)
      h5writeAttribute(as.integer(values), gid, attr_name)
      H5Gclose(gid)
      H5Fclose(fid)
    }, error = function(e3) {
      warning("Could not write shape attribute via rhdf5: ", e3$message)
    })
    
    # Use Python/h5py to fix the attribute type if reticulate available
    if (requireNamespace("reticulate", quietly = TRUE)) {
      fix_shape_with_python(h5file, group_path, attr_name, values)
    }
    
    return(invisible(TRUE))
  })
}


#' Fix shape attribute using h5py (ensures proper int64 numpy array)
fix_shape_with_python <- function(h5file, group_path, attr_name, values) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("reticulate not available - shape attribute may not be readable by anndata")
    return(invisible(FALSE))
  }
  
  tryCatch({
    py_code <- sprintf('
import h5py
import numpy as np

with h5py.File("%s", "a") as f:
    grp = f["%s"]
    # Delete existing attribute if present
    if "%s" in grp.attrs:
        del grp.attrs["%s"]
    # Write as proper int64 numpy array
    grp.attrs["%s"] = np.array([%s], dtype=np.int64)
', h5file, group_path, attr_name, attr_name, attr_name, paste(values, collapse = ", "))
    
    reticulate::py_run_string(py_code)
    return(invisible(TRUE))
    
  }, error = function(e) {
    warning("Could not fix shape with Python: ", e$message)
    return(invisible(FALSE))
  })
}


#' Finalize H5AD file using Python/h5py to ensure correct encoding
#' This fixes all encoding attributes to be compatible with anndata
finalize_h5ad_with_python <- function(h5file, n_cells, n_genes) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("reticulate not available - H5AD may not be readable by anndata")
    return(invisible(FALSE))
  }
  
  log_message("  Finalizing H5AD encoding with Python...")
  
  tryCatch({
    py_code <- sprintf('
import h5py
import numpy as np

with h5py.File("%s", "a") as f:
    # Fix X group encoding
    if "X" in f:
        x = f["X"]
        for attr in list(x.attrs.keys()):
            del x.attrs[attr]
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["encoding-version"] = "0.1.0"
        x.attrs.create("shape", data=np.array([%d, %d], dtype=np.int64))
    
    # Fix obs group and all its columns
    if "obs" in f:
        obs = f["obs"]
        
        # Get column-order if exists
        col_order = None
        if "column-order" in obs.attrs:
            col_order_raw = obs.attrs["column-order"]
            if hasattr(col_order_raw, "tolist"):
                col_order = [c.decode("utf-8") if isinstance(c, bytes) else c for c in col_order_raw.tolist()]
            elif isinstance(col_order_raw, (list, tuple)):
                col_order = [c.decode("utf-8") if isinstance(c, bytes) else c for c in col_order_raw]
        
        # Clear all attributes on obs group
        for attr in list(obs.attrs.keys()):
            del obs.attrs[attr]
        
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "_index"
        
        if col_order and len(col_order) > 0:
            obs.attrs.create("column-order", data=np.array(col_order, dtype="S"))
        
        # Fix encoding on each column in obs
        for col_name in obs.keys():
            if col_name == "_index":
                continue
            col = obs[col_name]
            
            if isinstance(col, h5py.Group):
                for attr in list(col.attrs.keys()):
                    del col.attrs[attr]
                col.attrs["encoding-type"] = "categorical"
                col.attrs["encoding-version"] = "0.2.0"
                col.attrs["ordered"] = False
            elif isinstance(col, h5py.Dataset):
                for attr in list(col.attrs.keys()):
                    del col.attrs[attr]
                if col.dtype.kind in ["U", "S", "O"]:
                    col.attrs["encoding-type"] = "string-array"
                    col.attrs["encoding-version"] = "0.2.0"
                else:
                    col.attrs["encoding-type"] = "array"
                    col.attrs["encoding-version"] = "0.2.0"
    
    # Fix var group and its columns
    if "var" in f:
        var = f["var"]
        
        # Get column order from keys (excluding _index)
        var_cols = [k for k in var.keys() if k != "_index"]
        
        for attr in list(var.attrs.keys()):
            del var.attrs[attr]
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "_index"
        
        # CRITICAL: Always write column-order, even if empty
        if var_cols:
            var.attrs.create("column-order", data=np.array(var_cols, dtype="S"))
        else:
            var.attrs.create("column-order", data=np.array([], dtype="S"))
        
        for col_name in var.keys():
            if col_name == "_index":
                continue
            col = var[col_name]
            if isinstance(col, h5py.Group):
                for attr in list(col.attrs.keys()):
                    del col.attrs[attr]
                col.attrs["encoding-type"] = "categorical"
                col.attrs["encoding-version"] = "0.2.0"
                col.attrs["ordered"] = False
            elif isinstance(col, h5py.Dataset):
                for attr in list(col.attrs.keys()):
                    del col.attrs[attr]
                if col.dtype.kind in ["U", "S", "O"]:
                    col.attrs["encoding-type"] = "string-array"
                    col.attrs["encoding-version"] = "0.2.0"
                else:
                    col.attrs["encoding-type"] = "array"
                    col.attrs["encoding-version"] = "0.2.0"
    
    # Fix empty dict groups and obsm contents
    for grp_name in ["obsm", "varm", "obsp", "varp", "uns", "layers"]:
        if grp_name in f:
            grp = f[grp_name]
            for attr in list(grp.attrs.keys()):
                del grp.attrs[attr]
            grp.attrs["encoding-type"] = "dict"
            grp.attrs["encoding-version"] = "0.1.0"
    
    # Fix obsm entries - check dimensions and remove bad ones
    if "obsm" in f:
        obsm = f["obsm"]
        n_obs = %d  # n_cells passed to function
        for key in list(obsm.keys()):
            item = obsm[key]
            if isinstance(item, h5py.Dataset):
                shape = item.shape
                if len(shape) == 1:
                    print(f"  Warning: obsm/{key} is 1D with shape {shape}, deleting")
                    del obsm[key]
                elif shape[0] != n_obs:
                    print(f"  Warning: obsm/{key} has shape {shape}, expected ({n_obs}, ...), deleting")
                    del obsm[key]

print("  H5AD encoding finalized successfully")
', h5file, n_cells, n_genes, n_cells)
    
    reticulate::py_run_string(py_code)
    return(invisible(TRUE))
    
  }, error = function(e) {
    warning("Could not finalize H5AD with Python: ", e$message)
    return(invisible(FALSE))
  })
}


# =============================================================================
# MAIN EXPORT FUNCTION
# =============================================================================

#' Main export function with chunked support for massive matrices
#' 
#' @param config_path Path to config.yaml
#' @param chunk_size Number of cells per chunk (default: 50000)
#' @param create_h5mu Whether to create H5MU with AIRR data (default: TRUE)
#' @param use_reticulate Use reticulate to create H5MU (default: TRUE)
#' @param verify Run verification after export (default: FALSE)
export_to_scanpy <- function(config_path = "config.yaml",
                             chunk_size = 50000,
                             create_h5mu = TRUE,
                             use_reticulate = TRUE,
                             verify = FALSE) {
  
  log_message("=== Starting Dual-Format Export (R + Python) ===")
  
  config <- yaml::read_yaml(config_path)
  input_path <- file.path(config$paths$results_dir, "FINAL_integrated_object.rds")
  export_dir <- file.path(config$paths$results_dir, "09_scanpy_export")
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!file.exists(input_path)) stop("Final object not found: ", input_path)
  
  log_message("Loading Seurat object...")
  obj <- readRDS(input_path)
  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  log_message(sprintf("Loaded object: %s cells x %s genes", 
                      format(n_cells, big.mark = ","),
                      format(n_genes, big.mark = ",")))
  
  # Get counts matrix
  counts <- obj@assays$RNA@layers$counts
  if (is.null(counts)) counts <- obj@assays$RNA@counts
  
  # Estimate nnz
  nnz_estimate <- tryCatch({
    sum(BPCells::matrix_stats(counts)$col_stats["nonzero",])
  }, error = function(e) {
    n_cells * 2000
  })
  
  log_message(sprintf("Estimated non-zeros: %s", format(nnz_estimate, big.mark = ",")))
  
  # Check if we need chunked export
  needs_chunking <- nnz_estimate > 2e9
  
  if (needs_chunking) {
    log_message("Matrix exceeds dgCMatrix limit - using chunked HDF5 export")
    h5ad_path <- export_chunked_h5ad(obj, export_dir, chunk_size)
  } else {
    log_message("Matrix within limits - using standard export")
    h5ad_path <- export_standard_h5ad(obj, export_dir)
  }
  
  # Export AIRR data
  airr_path <- export_airr_data(obj, export_dir)
  
  # Create H5MU
  h5mu_path <- NULL
  if (create_h5mu && !is.null(airr_path)) {
    log_message("\n[3/4] Creating H5MU (MuData) with scirpy-compatible AIRR...")
    
    if (use_reticulate && requireNamespace("reticulate", quietly = TRUE)) {
      h5mu_path <- create_h5mu_with_reticulate(h5ad_path, airr_path, export_dir)
    } else {
      generate_h5mu_creator_script(export_dir)
      log_message("  Generated create_h5mu.py - run manually or with: python3 create_h5mu.py")
    }
  } else {
    log_message("\n[3/4] Skipping H5MU creation (no AIRR data or create_h5mu=FALSE)")
  }
  
  log_message("\n[4/4] Generating documentation and verification scripts...")
  generate_verification_script(export_dir)
  generate_h5mu_creator_script(export_dir)
  
  log_message("=== Export Complete ===")
  log_message(sprintf("Output directory: %s", export_dir))
  
  return(invisible(export_dir))
}


# =============================================================================
# H5AD EXPORT FUNCTIONS
# =============================================================================

#' Export large matrix in chunks directly to H5AD format
export_chunked_h5ad <- function(obj, export_dir, chunk_size = 50000) {
  
  h5ad_path <- file.path(export_dir, "adata.h5ad")
  
  if (file.exists(h5ad_path)) unlink(h5ad_path)
  
  log_message("Creating H5AD file with chunked matrix export...")
  
  counts <- obj@assays$RNA@layers$counts
  if (is.null(counts)) counts <- obj@assays$RNA@counts
  
  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  cell_names <- colnames(obj)
  gene_names <- rownames(obj)
  
  # Create HDF5 file
  h5createFile(h5ad_path)
  
  # Create all required groups
  h5createGroup(h5ad_path, "X")
  h5createGroup(h5ad_path, "obs")
  h5createGroup(h5ad_path, "var")
  h5createGroup(h5ad_path, "obsm")
  h5createGroup(h5ad_path, "varm")
  h5createGroup(h5ad_path, "obsp")
  h5createGroup(h5ad_path, "varp")
  h5createGroup(h5ad_path, "uns")
  h5createGroup(h5ad_path, "layers")
  
  # Process matrix in chunks
  n_chunks <- ceiling(n_cells / chunk_size)
  log_message(sprintf("Processing %d chunks of ~%d cells each", n_chunks, chunk_size))
  
  all_data <- list()
  all_indices <- list()
  indptr_vec <- 0L
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_cells)
    
    if (i %% 5 == 1 || i == n_chunks) {
      log_message(sprintf("  Chunk %d/%d: cells %d-%d", i, n_chunks, start_idx, end_idx))
    }
    
    chunk_counts <- counts[, start_idx:end_idx, drop = FALSE]
    chunk_sparse <- as(chunk_counts, "dgCMatrix")
    
    # Transpose to cells x genes (CSR format)
    chunk_csr <- as(t(chunk_sparse), "RsparseMatrix")
    
    all_data[[i]] <- chunk_csr@x
    all_indices[[i]] <- chunk_csr@j
    
    # Update indptr
    chunk_indptr <- chunk_csr@p[-1]
    current_offset <- indptr_vec[length(indptr_vec)]
    indptr_vec <- c(indptr_vec, chunk_indptr + current_offset)
    
    rm(chunk_sparse, chunk_csr, chunk_counts)
    if (i %% 10 == 0) gc()
  }
  
  log_message("Combining chunks and writing to H5AD...")
  
  combined_data <- unlist(all_data)
  combined_indices <- unlist(all_indices)
  rm(all_data, all_indices); gc()
  
  # Write sparse matrix components
  h5write(combined_data, h5ad_path, "X/data")
  h5write(as.integer(combined_indices), h5ad_path, "X/indices")
  h5write(as.integer(indptr_vec), h5ad_path, "X/indptr")
  
  rm(combined_data, combined_indices, indptr_vec); gc()
  
  # Write metadata
  log_message("Writing cell metadata...")
  write_obs_to_h5ad(obj@meta.data, h5ad_path, cell_names)
  
  log_message("Writing gene metadata...")
  write_var_to_h5ad(gene_names, h5ad_path)
  
  log_message("Writing embeddings...")
  write_obsm_to_h5ad(obj, h5ad_path)
  
  # CRITICAL: Use Python to finalize encoding attributes
  # This ensures the H5AD is readable by anndata/scanpy
  if (requireNamespace("reticulate", quietly = TRUE)) {
    finalize_h5ad_with_python(h5ad_path, n_cells, n_genes)
  } else {
    # Fallback: try to set attributes with rhdf5
    write_int_array_attribute(h5ad_path, "X", "shape", c(n_cells, n_genes))
    fid <- H5Fopen(h5ad_path)
    gid <- H5Gopen(fid, "X")
    h5writeAttribute("csr_matrix", gid, "encoding-type")
    h5writeAttribute("0.1.0", gid, "encoding-version")
    H5Gclose(gid)
    H5Fclose(fid)
    write_empty_groups(h5ad_path)
  }
  
  log_message(sprintf("H5AD saved: %s", h5ad_path))
  return(h5ad_path)
}


#' Standard export for smaller matrices (< 2 billion non-zeros)
export_standard_h5ad <- function(obj, export_dir) {
  
  h5ad_path <- file.path(export_dir, "adata.h5ad")
  if (file.exists(h5ad_path)) unlink(h5ad_path)
  
  log_message("Creating H5AD file...")
  
  counts <- obj@assays$RNA@layers$counts
  if (is.null(counts)) counts <- obj@assays$RNA@counts
  
  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  cell_names <- colnames(obj)
  gene_names <- rownames(obj)
  
  # Convert to sparse CSR
  mat_sparse <- as(counts, "dgCMatrix")
  mat_csr <- as(t(mat_sparse), "RsparseMatrix")
  
  # Create file and groups
  h5createFile(h5ad_path)
  h5createGroup(h5ad_path, "X")
  h5createGroup(h5ad_path, "obs")
  h5createGroup(h5ad_path, "var")
  h5createGroup(h5ad_path, "obsm")
  h5createGroup(h5ad_path, "varm")
  h5createGroup(h5ad_path, "obsp")
  h5createGroup(h5ad_path, "varp")
  h5createGroup(h5ad_path, "uns")
  h5createGroup(h5ad_path, "layers")
  
  # Write sparse matrix
  h5write(mat_csr@x, h5ad_path, "X/data")
  h5write(as.integer(mat_csr@j), h5ad_path, "X/indices")
  h5write(as.integer(mat_csr@p), h5ad_path, "X/indptr")
  
  # Write metadata
  write_obs_to_h5ad(obj@meta.data, h5ad_path, cell_names)
  write_var_to_h5ad(gene_names, h5ad_path)
  write_obsm_to_h5ad(obj, h5ad_path)
  
  # CRITICAL: Use Python to finalize encoding attributes
  if (requireNamespace("reticulate", quietly = TRUE)) {
    finalize_h5ad_with_python(h5ad_path, n_cells, n_genes)
  } else {
    write_int_array_attribute(h5ad_path, "X", "shape", c(n_cells, n_genes))
    fid <- H5Fopen(h5ad_path)
    gid <- H5Gopen(fid, "X")
    h5writeAttribute("csr_matrix", gid, "encoding-type")
    h5writeAttribute("0.1.0", gid, "encoding-version")
    H5Gclose(gid)
    H5Fclose(fid)
    write_empty_groups(h5ad_path)
  }
  
  log_message(sprintf("H5AD saved: %s", h5ad_path))
  return(h5ad_path)
}


# =============================================================================
# H5AD COMPONENT WRITERS
# =============================================================================

#' Write obs (cell metadata) to H5AD
write_obs_to_h5ad <- function(meta_data, h5ad_path, cell_names) {
  
  # Write cell barcodes as index
  h5write(cell_names, h5ad_path, "obs/_index")
  
  # Set attributes properly
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "obs")
  
  h5writeAttribute("_index", gid, "_index")
  
  # FIX: Only write column-order if there are columns
  # Empty character(0) causes HDF5 "Bad value" error
  col_names <- colnames(meta_data)
  if (length(col_names) > 0) {
    h5writeAttribute(col_names, gid, "column-order")
  }
  # If empty, simply omit - AnnData handles missing column-order
  
  h5writeAttribute("dataframe", gid, "encoding-type")
  h5writeAttribute("0.2.0", gid, "encoding-version")
  
  H5Gclose(gid)
  H5Fclose(fid)
  
  # Write each metadata column
  for (col in colnames(meta_data)) {
    values <- meta_data[[col]]
    col_path <- paste0("obs/", col)
    
    tryCatch({
      if (is.factor(values)) {
        # Categorical encoding
        h5createGroup(h5ad_path, col_path)
        
        codes <- as.integer(values) - 1L  # 0-indexed
        codes[is.na(values)] <- -1L  # NA handling
        
        h5write(codes, h5ad_path, paste0(col_path, "/codes"))
        h5write(levels(values), h5ad_path, paste0(col_path, "/categories"))
        
        fid <- H5Fopen(h5ad_path)
        gid <- H5Gopen(fid, col_path)
        h5writeAttribute("categorical", gid, "encoding-type")
        h5writeAttribute("0.2.0", gid, "encoding-version")
        # Use integer 0 instead of FALSE to avoid potential issues
        h5writeAttribute(0L, gid, "ordered")
        H5Gclose(gid)
        H5Fclose(fid)
        
      } else if (is.character(values)) {
        # String array - handle NAs
        values[is.na(values)] <- ""
        h5write(values, h5ad_path, col_path)
        
        # Add string encoding attribute
        fid <- H5Fopen(h5ad_path)
        did <- H5Dopen(fid, col_path)
        h5writeAttribute("string-array", did, "encoding-type")
        h5writeAttribute("0.2.0", did, "encoding-version")
        H5Dclose(did)
        H5Fclose(fid)
        
      } else if (is.numeric(values)) {
        # Numeric array
        h5write(as.numeric(values), h5ad_path, col_path)
        
        # Add array encoding attribute
        fid <- H5Fopen(h5ad_path)
        did <- H5Dopen(fid, col_path)
        h5writeAttribute("array", did, "encoding-type")
        h5writeAttribute("0.2.0", did, "encoding-version")
        H5Dclose(did)
        H5Fclose(fid)
        
      } else if (is.logical(values)) {
        # Boolean as integer
        h5write(as.integer(values), h5ad_path, col_path)
        
        fid <- H5Fopen(h5ad_path)
        did <- H5Dopen(fid, col_path)
        h5writeAttribute("array", did, "encoding-type")
        h5writeAttribute("0.2.0", did, "encoding-version")
        H5Dclose(did)
        H5Fclose(fid)
      }
    }, error = function(e) {
      warning(sprintf("Could not write column '%s': %s", col, e$message))
    })
  }
}


#' Write var (gene metadata) to H5AD
write_var_to_h5ad <- function(gene_names, h5ad_path) {
  
  h5write(gene_names, h5ad_path, "var/_index")
  
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "var")
  
  h5writeAttribute("_index", gid, "_index")
  # NOTE: Omit column-order when empty - AnnData handles missing attribute fine
  # Empty character(0) causes HDF5 error
  h5writeAttribute("dataframe", gid, "encoding-type")
  h5writeAttribute("0.2.0", gid, "encoding-version")
  
  H5Gclose(gid)
  H5Fclose(fid)
}


#' Write obsm (embeddings) to H5AD
write_obsm_to_h5ad <- function(obj, h5ad_path) {
  
  n_cells <- ncol(obj)
  
  for (red_name in names(obj@reductions)) {
    emb <- Embeddings(obj, reduction = red_name)
    
    # Ensure it's a matrix with correct dimensions
    if (!is.matrix(emb)) {
      emb <- as.matrix(emb)
    }
    
    # Check dimensions - should be n_cells x n_dims
    if (nrow(emb) != n_cells) {
      warning(sprintf("Embedding %s has %d rows but expected %d cells. Skipping.", 
                      red_name, nrow(emb), n_cells))
      next
    }
    
    # AnnData convention: X_pca, X_umap, etc.
    key <- paste0("X_", tolower(red_name))
    
    log_message(sprintf("  Writing %s: %d x %d", key, nrow(emb), ncol(emb)))
    
    # Write as 2D array - force dimensions
    h5write(emb, h5ad_path, paste0("obsm/", key))
  }
  
  # Set encoding for obsm group
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "obsm")
  h5writeAttribute("dict", gid, "encoding-type")
  h5writeAttribute("0.1.0", gid, "encoding-version")
  H5Gclose(gid)
  H5Fclose(fid)
}


#' Write empty groups with proper encoding
write_empty_groups <- function(h5ad_path) {
  
  empty_dict_groups <- c("varm", "obsp", "varp", "uns", "layers")
  
  for (grp in empty_dict_groups) {
    fid <- H5Fopen(h5ad_path)
    gid <- H5Gopen(fid, grp)
    h5writeAttribute("dict", gid, "encoding-type")
    h5writeAttribute("0.1.0", gid, "encoding-version")
    H5Gclose(gid)
    H5Fclose(fid)
  }
}


# =============================================================================
# AIRR/VDJ EXPORT
# =============================================================================

#' Export AIRR data from scRepertoire columns
export_airr_data <- function(obj, export_dir) {
  
  has_vdj <- any(c("CTaa", "CTgene", "TCR1", "IGH") %in% colnames(obj@meta.data))
  
  if (!has_vdj) {
    log_message("  No VDJ data found in metadata")
    return(NULL)
  }
  
  log_message("Exporting VDJ data in AIRR format...")
  airr_df <- parse_screpertoire_to_airr(obj@meta.data)
  
  if (!is.null(airr_df) && nrow(airr_df) > 0) {
    airr_path <- file.path(export_dir, "airr_rearrangement.tsv")
    fwrite(airr_df, file = airr_path, sep = "\t", row.names = FALSE)
    log_message(sprintf("  Exported %s AIRR chains from %s cells", 
                        format(nrow(airr_df), big.mark = ","),
                        format(length(unique(airr_df$cell_id)), big.mark = ",")))
    return(airr_path)
  }
  
  log_message("  Could not parse VDJ data into AIRR format")
  return(NULL)
}


#' Parse scRepertoire columns to AIRR format
parse_screpertoire_to_airr <- function(meta_data) {
  
  has_individual_chains <- any(c("TCR1", "TCR2", "IGH", "IGLC") %in% colnames(meta_data))
  has_combined <- "CTgene" %in% colnames(meta_data)
  
  if (!has_individual_chains && !has_combined) return(NULL)
  
  nested_rows <- lapply(seq_len(nrow(meta_data)), function(i) {
    cell_id <- rownames(meta_data)[i]
    row <- meta_data[i, , drop = FALSE]
    cell_airr_rows <- list()
    
    if (has_individual_chains) {
      chain_specs <- list(
        list(gene_col = "TCR1", aa_col = "cdr3_aa1", nt_col = "cdr3_nt1"),
        list(gene_col = "TCR2", aa_col = "cdr3_aa2", nt_col = "cdr3_nt2"),
        list(gene_col = "IGH", aa_col = "cdr3_aa1", nt_col = "cdr3_nt1"),
        list(gene_col = "IGLC", aa_col = "cdr3_aa2", nt_col = "cdr3_nt2")
      )
      
      for (spec in chain_specs) {
        gene_col <- spec$gene_col
        if (gene_col %in% colnames(meta_data) && !is.na(row[[gene_col]])) {
          gene_str <- as.character(row[[gene_col]])
          if (gene_str != "" && gene_str != "NA") {
            airr_row <- parse_gene_string_to_airr(
              cell_id = cell_id,
              gene_str = gene_str,
              cdr3_aa = if(spec$aa_col %in% colnames(meta_data)) as.character(row[[spec$aa_col]]) else NA,
              cdr3_nt = if(spec$nt_col %in% colnames(meta_data)) as.character(row[[spec$nt_col]]) else NA
            )
            if (!is.null(airr_row)) cell_airr_rows[[length(cell_airr_rows) + 1]] <- airr_row
          }
        }
      }
    }
    
    if (length(cell_airr_rows) == 0 && has_combined && !is.na(row[["CTgene"]])) {
      ct_gene <- as.character(row[["CTgene"]])
      ct_aa <- if("CTaa" %in% colnames(meta_data)) as.character(row[["CTaa"]]) else ""
      ct_nt <- if("CTnt" %in% colnames(meta_data)) as.character(row[["CTnt"]]) else ""
      
      if (!is.na(ct_gene) && ct_gene != "" && ct_gene != "NA") {
        gene_parts <- strsplit(ct_gene, "_")[[1]]
        aa_parts <- if(!is.na(ct_aa) && ct_aa != "") strsplit(ct_aa, "_")[[1]] else rep("", length(gene_parts))
        nt_parts <- if(!is.na(ct_nt) && ct_nt != "") strsplit(ct_nt, "_")[[1]] else rep("", length(gene_parts))
        
        for (j in seq_along(gene_parts)) {
          if (gene_parts[j] != "" && gene_parts[j] != "NA") {
            airr_row <- parse_gene_string_to_airr(
              cell_id = cell_id,
              gene_str = gene_parts[j],
              cdr3_aa = if(j <= length(aa_parts)) aa_parts[j] else NA,
              cdr3_nt = if(j <= length(nt_parts)) nt_parts[j] else NA
            )
            if (!is.null(airr_row)) cell_airr_rows[[length(cell_airr_rows) + 1]] <- airr_row
          }
        }
      }
    }
    
    return(cell_airr_rows)
  })
  
  airr_rows <- unlist(nested_rows, recursive = FALSE)
  
  if (!is.null(airr_rows) && length(airr_rows) > 0) {
    return(do.call(rbind, lapply(airr_rows, as.data.frame)))
  }
  return(NULL)
}


#' Parse gene string to AIRR row
parse_gene_string_to_airr <- function(cell_id, gene_str, cdr3_aa = NA, cdr3_nt = NA) {
  genes <- strsplit(gene_str, "\\.")[[1]]
  
  if (length(genes) == 0 || all(genes == "None") || all(genes == "NA")) return(NULL)
  
  airr_row <- list(
    cell_id = cell_id,
    sequence_id = paste0(cell_id, "_", gsub("\\.", "-", gene_str)),
    locus = NA_character_,
    productive = TRUE,
    v_call = NA_character_,
    d_call = NA_character_,
    j_call = NA_character_,
    c_call = NA_character_,
    junction_aa = NA_character_,
    junction = NA_character_
  )
  
  v_gene <- genes[1]
  if (!is.na(v_gene) && v_gene != "None" && v_gene != "") {
    airr_row$v_call <- v_gene
    
    if (grepl("^TRA", v_gene)) airr_row$locus <- "TRA"
    else if (grepl("^TRB", v_gene)) airr_row$locus <- "TRB"
    else if (grepl("^TRG", v_gene)) airr_row$locus <- "TRG"
    else if (grepl("^TRD", v_gene)) airr_row$locus <- "TRD"
    else if (grepl("^IGH", v_gene)) airr_row$locus <- "IGH"
    else if (grepl("^IGK", v_gene)) airr_row$locus <- "IGK"
    else if (grepl("^IGL", v_gene)) airr_row$locus <- "IGL"
  }
  
  for (gene in genes[-1]) {
    if (!is.na(gene) && gene != "None" && gene != "") {
      if (grepl("^TR[ABDG]J|^IG[HKL]J", gene)) airr_row$j_call <- gene
      else if (grepl("^TR[BD]D|^IGHD", gene)) airr_row$d_call <- gene
      else if (grepl("^TR[ABDG]C|^IG[HKL]C|^IGH[MGDEA]", gene)) airr_row$c_call <- gene
    }
  }
  
  if (!is.na(cdr3_aa) && cdr3_aa != "" && cdr3_aa != "NA") airr_row$junction_aa <- cdr3_aa
  if (!is.na(cdr3_nt) && cdr3_nt != "" && cdr3_nt != "NA") airr_row$junction <- cdr3_nt
  
  if (!is.na(airr_row$locus)) return(airr_row)
  return(NULL)
}


# =============================================================================
# H5MU CREATION
# =============================================================================

#' Create H5MU using reticulate
create_h5mu_with_reticulate <- function(h5ad_path, airr_path, export_dir) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    log_message("  reticulate not available")
    return(NULL)
  }
  
  library(reticulate)
  
  # Check Python packages
  required_pkgs <- c("scanpy", "anndata", "muon", "scirpy")
  for (pkg in required_pkgs) {
    if (!py_module_available(pkg)) {
      log_message(sprintf("  Missing Python package: %s", pkg))
      return(NULL)
    }
  }
  
  h5mu_path <- file.path(export_dir, "mudata.h5mu")
  
  tryCatch({
    log_message("  Creating MuData via reticulate...")
    
    py_code <- sprintf('
import scanpy as sc
import scirpy as ir
import muon as mu

# Load GEX
adata_gex = sc.read_h5ad("%s")
print(f"  GEX: {adata_gex.n_obs} cells x {adata_gex.n_vars} genes")

# Load AIRR
adata_airr = ir.io.read_airr("%s")
print(f"  AIRR: {adata_airr.n_obs} cells")

# Create MuData
mdata = mu.MuData({"gex": adata_gex, "airr": adata_airr})

# Run scirpy preprocessing
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

# Save
mdata.write("%s")
print(f"  Saved: %s")
', h5ad_path, airr_path, h5mu_path, h5mu_path)
    
    py_run_string(py_code)
    
    log_message("  H5MU created successfully with scirpy chain indexing")
    return(h5mu_path)
    
  }, error = function(e) {
    log_message(sprintf("  Error creating H5MU: %s", e$message))
    log_message("  Generating Python script for manual creation...")
    generate_h5mu_creator_script(export_dir)
    return(NULL)
  })
}


# =============================================================================
# SCRIPT GENERATORS
# =============================================================================

#' Generate Python script to create H5MU
generate_h5mu_creator_script <- function(export_dir) {
  
  script <- '#!/usr/bin/env python3
"""
Create H5MU (MuData) file combining GEX and AIRR data.
Run this if automatic creation via reticulate failed.

Usage:
    python3 create_h5mu.py
"""

import sys
try:
    import scanpy as sc
    import scirpy as ir
    import muon as mu
except ImportError as e:
    print(f"Missing package: {e}")
    print("Install with: pip install scanpy scirpy muon")
    sys.exit(1)

print("Creating MuData (H5MU) file...")

# Load GEX
print("Loading gene expression data...")
adata_gex = sc.read_h5ad("adata.h5ad")
print(f"  GEX: {adata_gex.n_obs} cells x {adata_gex.n_vars} genes")

# Load AIRR
print("Loading AIRR data...")
adata_airr = ir.io.read_airr("airr_rearrangement.tsv")
print(f"  AIRR: {adata_airr.n_obs} cells with receptor data")

# Create MuData
print("Creating MuData object...")
mdata = mu.MuData({"gex": adata_gex, "airr": adata_airr})

# Run scirpy preprocessing
print("Running scirpy chain QC...")
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

# Show summary
if "receptor_type" in mdata["airr"].obs.columns:
    print("\\nReceptor type distribution:")
    print(mdata["airr"].obs["receptor_type"].value_counts())

# Save
mdata.write("mudata.h5mu")
print("\\nSaved: mudata.h5mu")
'
  
  writeLines(script, file.path(export_dir, "create_h5mu.py"))
  log_message("  Generated: create_h5mu.py")
}


#' Generate verification script
generate_verification_script <- function(export_dir) {
  
  script <- '#!/usr/bin/env python3
"""
Verify the exported H5AD file loads correctly in scanpy.
"""

import sys

try:
    import scanpy as sc
    import anndata as ad
except ImportError as e:
    print(f"Missing: {e}")
    sys.exit(1)

print(f"scanpy: {sc.__version__}")
print(f"anndata: {ad.__version__}")

print("\\nLoading adata.h5ad...")
try:
    adata = sc.read_h5ad("adata.h5ad")
    print(f"✓ Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  obs columns: {list(adata.obs.columns)[:5]}...")
    print(f"  obsm keys: {list(adata.obsm.keys())}")
    
    if hasattr(adata.X, "nnz"):
        print(f"  Matrix: sparse, {adata.X.nnz:,} non-zeros")
    
    print("\\n✓ H5AD verification passed!")
    
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'
  
  writeLines(script, file.path(export_dir, "verify_export.py"))
  log_message("  Generated: verify_export.py")
}
