# R/80_export_to_scanpy.R
# 
# Unified export for both R and Python users
# Handles massive matrices (>2^31 non-zeros) via chunked export

suppressPackageStartupMessages({
  
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(yaml)
  library(BPCells)
  library(rhdf5)  # For direct HDF5 writing
})

source("R/00_utils.R")

#' Main export function with chunked support for massive matrices
#' 
#' @param config_path Path to config.yaml
#' @param chunk_size Number of cells per chunk (default: 50000)
export_to_scanpy <- function(config_path = "config.yaml",
                             chunk_size = 50000) {
  
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
    # BPCells can give us this info
    sum(BPCells::matrix_stats(counts)$col_stats["nonzero",])
  }, error = function(e) {
    n_cells * 2000  # Rough estimate: ~2000 genes per cell
  })
  
  log_message(sprintf("Estimated non-zeros: %s", format(nnz_estimate, big.mark = ",")))
  
  # Check if we need chunked export
  needs_chunking <- nnz_estimate > 2e9  # 2 billion threshold
  
  if (needs_chunking) {
    log_message("Matrix exceeds dgCMatrix limit - using chunked HDF5 export")
    h5ad_path <- export_chunked_h5ad(obj, export_dir, chunk_size)
  } else {
    log_message("Matrix within limits - using standard export")
    h5ad_path <- export_standard_h5ad(obj, export_dir)
  }
  
  # Export AIRR data
  airr_path <- export_airr_data(obj, export_dir)
  
  # Generate Python scripts
  generate_helper_scripts(export_dir, has_airr = !is.null(airr_path))
  
  # Generate README
  generate_readme(export_dir, h5ad_path, airr_path, chunked = needs_chunking)
  
  log_message("=== Export Complete ===")
  log_message(sprintf("Output directory: %s", export_dir))
  
  return(invisible(export_dir))
}


#' Export massive matrix in chunks directly to H5AD format
#' 
#' This writes directly to HDF5 without creating a dgCMatrix,
#' bypassing the 2^31 non-zero limit.
export_chunked_h5ad <- function(obj, export_dir, chunk_size = 50000) {
  
  h5ad_path <- file.path(export_dir, "adata.h5ad")
  
  # Remove existing file
  if (file.exists(h5ad_path)) unlink(h5ad_path)
  
  log_message("Creating H5AD file with chunked matrix export...")
  
  # Get dimensions
  counts <- obj@assays$RNA@layers$counts
  if (is.null(counts)) counts <- obj@assays$RNA@counts
  
  n_cells <- ncol(obj)
  n_genes <- nrow(obj)
  cell_names <- colnames(obj)
  gene_names <- rownames(obj)
  
  # Create HDF5 file with H5AD structure
  h5createFile(h5ad_path)
  
  # Calculate chunks
  n_chunks <- ceiling(n_cells / chunk_size)
  log_message(sprintf("Processing %d chunks of ~%d cells each", n_chunks, chunk_size))
  
  # We'll collect CSR components across chunks
  all_data <- list()
  all_indices <- list()
  indptr <- 0L  # Running pointer
  indptr_vec <- c(0L)
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_cells)
    
    log_message(sprintf("  Chunk %d/%d: cells %d-%d", i, n_chunks, start_idx, end_idx))
    
    # Extract chunk using BPCells (memory-efficient)
    chunk_cells <- cell_names[start_idx:end_idx]
    chunk_counts <- counts[, start_idx:end_idx, drop = FALSE]
    
    # Convert this smaller chunk to dgCMatrix
    chunk_sparse <- as(chunk_counts, "dgCMatrix")
    
    # Extract CSR components (we need to transpose for AnnData which is cells x genes)
    # dgCMatrix is CSC (genes x cells), we need CSR (cells x genes)
    chunk_csr <- as(t(chunk_sparse), "RsparseMatrix")
    
    # Collect data
    all_data[[i]] <- chunk_csr@x
    all_indices[[i]] <- chunk_csr@j  # Column indices (genes)
    
    # Update indptr (row pointers for cells)
    chunk_indptr <- chunk_csr@p[-1]  # Remove first 0
    indptr_vec <- c(indptr_vec, chunk_indptr + indptr)
    indptr <- indptr_vec[length(indptr_vec)]
    
    rm(chunk_sparse, chunk_csr, chunk_counts)
    gc()
  }
  
  # Combine all chunks
  log_message("Combining chunks and writing to H5AD...")
  
  combined_data <- unlist(all_data)
  combined_indices <- unlist(all_indices)
  
  rm(all_data, all_indices)
  gc()
  
  # Write to HDF5 in AnnData format
  # AnnData stores X in CSR format under /X with data, indices, indptr
  
  # Create groups
  h5createGroup(h5ad_path, "X")
  h5createGroup(h5ad_path, "obs")
  h5createGroup(h5ad_path, "var")
  h5createGroup(h5ad_path, "obsm")
  h5createGroup(h5ad_path, "uns")
  
  # Write sparse matrix in CSR format
  h5write(combined_data, h5ad_path, "X/data")
  h5write(as.integer(combined_indices), h5ad_path, "X/indices")
  h5write(as.integer(indptr_vec), h5ad_path, "X/indptr")
  
  # Write shape attribute
  h5write(c(n_cells, n_genes), h5ad_path, "X/shape")
  
  # Add encoding metadata for AnnData
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "X")
  h5writeAttribute("csr_matrix", gid, "encoding-type")
  h5writeAttribute("0.1.0", gid, "encoding-version")
  H5Gclose(gid)
  H5Fclose(fid)
  
  rm(combined_data, combined_indices, indptr_vec)
  gc()
  
  # Write obs (cell metadata)
  log_message("Writing cell metadata...")
  write_obs_to_h5ad(obj@meta.data, h5ad_path, cell_names)
  
  # Write var (gene metadata)
  log_message("Writing gene metadata...")
  write_var_to_h5ad(gene_names, h5ad_path)
  
  # Write embeddings
  log_message("Writing embeddings...")
  write_obsm_to_h5ad(obj, h5ad_path)
  
  log_message(sprintf("H5AD saved: %s", h5ad_path))
  return(h5ad_path)
}


#' Write obs (cell metadata) to H5AD
write_obs_to_h5ad <- function(meta_data, h5ad_path, cell_names) {
  
  # Write index (cell barcodes)
  h5write(cell_names, h5ad_path, "obs/_index")
  
  # Add index attribute
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "obs")
  h5writeAttribute("_index", gid, "_index")
  h5writeAttribute(c("_index", colnames(meta_data)), gid, "column-order")
  h5writeAttribute("dataframe", gid, "encoding-type")
  h5writeAttribute("0.2.0", gid, "encoding-version")
  H5Gclose(gid)
  H5Fclose(fid)
  
  # Write each column
  for (col in colnames(meta_data)) {
    values <- meta_data[[col]]
    
    if (is.factor(values)) {
      # Categorical encoding for AnnData
      h5createGroup(h5ad_path, paste0("obs/", col))
      h5write(as.integer(values) - 1L, h5ad_path, paste0("obs/", col, "/codes"))
      h5write(levels(values), h5ad_path, paste0("obs/", col, "/categories"))
      
      # Add categorical attributes
      fid <- H5Fopen(h5ad_path)
      gid <- H5Gopen(fid, paste0("obs/", col))
      h5writeAttribute("categorical", gid, "encoding-type")
      h5writeAttribute("0.2.0", gid, "encoding-version")
      h5writeAttribute(FALSE, gid, "ordered")
      H5Gclose(gid)
      H5Fclose(fid)
      
    } else if (is.character(values)) {
      h5write(values, h5ad_path, paste0("obs/", col))
      
    } else if (is.numeric(values)) {
      h5write(as.numeric(values), h5ad_path, paste0("obs/", col))
      
    } else if (is.logical(values)) {
      h5write(as.integer(values), h5ad_path, paste0("obs/", col))
    }
  }
}


#' Write var (gene metadata) to H5AD
write_var_to_h5ad <- function(gene_names, h5ad_path) {
  
  # Write index (gene names)
  h5write(gene_names, h5ad_path, "var/_index")
  
  # Add attributes
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "var")
  h5writeAttribute("_index", gid, "_index")
  h5writeAttribute(c("_index"), gid, "column-order")
  h5writeAttribute("dataframe", gid, "encoding-type")
  h5writeAttribute("0.2.0", gid, "encoding-version")
  H5Gclose(gid)
  H5Fclose(fid)
}


#' Write obsm (embeddings) to H5AD
write_obsm_to_h5ad <- function(obj, h5ad_path) {
  
  for (red_name in names(obj@reductions)) {
    emb <- Embeddings(obj, reduction = red_name)
    key <- paste0("X_", red_name)
    
    h5write(emb, h5ad_path, paste0("obsm/", key))
    log_message(sprintf("  Wrote embedding: %s (%d dims)", key, ncol(emb)))
  }
  
  # Add obsm encoding
  fid <- H5Fopen(h5ad_path)
  gid <- H5Gopen(fid, "obsm")
  h5writeAttribute("dict", gid, "encoding-type")
  h5writeAttribute("0.1.0", gid, "encoding-version")
  H5Gclose(gid)
  H5Fclose(fid)
}


#' Standard export for smaller matrices
export_standard_h5ad <- function(obj, export_dir) {
  
  # Try available packages
  if (requireNamespace("anndataR", quietly = TRUE)) {
    library(anndataR)
    h5ad_path <- file.path(export_dir, "adata.h5ad")
    adata <- from_Seurat(obj)
    write_h5ad(adata, h5ad_path)
    return(h5ad_path)
  }
  
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    library(zellkonverter)
    library(SingleCellExperiment)
    h5ad_path <- file.path(export_dir, "adata.h5ad")
    sce <- as.SingleCellExperiment(obj)
    writeH5AD(sce, file = h5ad_path, X_name = "counts")
    return(h5ad_path)
  }
  
  # Fallback to our chunked method even for smaller matrices
  log_message("No standard H5AD package found, using direct HDF5 export")
  return(export_chunked_h5ad(obj, export_dir, chunk_size = 100000))
}


#' Export AIRR data from scRepertoire columns
export_airr_data <- function(obj, export_dir) {
  
  has_vdj <- any(c("CTaa", "CTgene", "TCR1", "IGH") %in% colnames(obj@meta.data))
  
  if (!has_vdj) {
    log_message("No VDJ data found in metadata")
    return(NULL)
  }
  
  log_message("Exporting VDJ data in AIRR format...")
  airr_df <- parse_screpertoire_to_airr(obj@meta.data)
  
  if (!is.null(airr_df) && nrow(airr_df) > 0) {
    airr_path <- file.path(export_dir, "airr_rearrangement.tsv")
    fwrite(airr_df, file = airr_path, sep = "\t", row.names = FALSE)
    log_message(sprintf("Exported %s AIRR chain records", 
                        format(nrow(airr_df), big.mark = ",")))
    return(airr_path)
  }
  
  return(NULL)
}


#' Generate helper Python scripts
generate_helper_scripts <- function(export_dir, has_airr = FALSE) {
  
  # Script to verify and optionally fix the H5AD
  verify_script <- '#!/usr/bin/env python3
"""
Verify and load the exported H5AD file.
Optionally creates MuData with AIRR data.

Usage:
    python3 load_data.py
"""

import os
import sys

try:
    import scanpy as sc
    import anndata as ad
    import pandas as pd
    import numpy as np
except ImportError as e:
    print(f"Missing package: {e}")
    print("Install with: pip install scanpy anndata pandas")
    sys.exit(1)

print(f"Scanpy version: {sc.__version__}")
print(f"AnnData version: {ad.__version__}")

# Load H5AD
print("\\nLoading adata.h5ad...")
try:
    adata = sc.read_h5ad("adata.h5ad")
    print(f"Successfully loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"\\nObservations (obs): {list(adata.obs.columns)[:10]}...")
    print(f"Embeddings (obsm): {list(adata.obsm.keys())}")
    
    # Basic QC
    print(f"\\nMatrix stats:")
    print(f"  Shape: {adata.X.shape}")
    print(f"  Non-zeros: {adata.X.nnz:,}")
    print(f"  Density: {adata.X.nnz / (adata.n_obs * adata.n_vars) * 100:.2f}%")
    
except Exception as e:
    print(f"Error loading H5AD: {e}")
    print("\\nTrying to diagnose...")
    
    import h5py
    with h5py.File("adata.h5ad", "r") as f:
        print("HDF5 structure:")
        def print_structure(name, obj):
            print(f"  {name}: {type(obj).__name__}")
        f.visititems(print_structure)
    sys.exit(1)

# Load AIRR and create MuData if available
if os.path.exists("airr_rearrangement.tsv"):
    print("\\n--- Loading AIRR data ---")
    try:
        import scirpy as ir
        import muon as mu
        
        adata_airr = ir.io.read_airr("airr_rearrangement.tsv")
        print(f"AIRR data: {adata_airr.n_obs} cells")
        
        # Create MuData
        print("\\nCreating MuData...")
        mdata = mu.MuData({"gex": adata, "airr": adata_airr})
        
        # Run scirpy QC
        ir.pp.index_chains(mdata)
        ir.tl.chain_qc(mdata)
        
        if "receptor_type" in mdata["airr"].obs.columns:
            print("\\nReceptor types:")
            print(mdata["airr"].obs["receptor_type"].value_counts())
        
        # Save MuData
        mdata.write("mudata.h5mu")
        print("\\nSaved: mudata.h5mu")
        
    except ImportError:
        print("scirpy/muon not installed - skipping MuData creation")
        print("Install with: pip install scirpy muon")
    except Exception as e:
        print(f"Error with AIRR data: {e}")

print("\\nDone!")
'
  
  writeLines(verify_script, file.path(export_dir, "load_data.py"))
  log_message("Generated load_data.py")
}


#' Generate README
generate_readme <- function(export_dir, h5ad_path, airr_path, chunked = FALSE) {
  
  readme <- sprintf('# Single-Cell Data Export

## Files

- `adata.h5ad` - AnnData object (gene expression, metadata, embeddings)
%s
- `load_data.py` - Python script to load and verify data

## Quick Start (Python)

```python
import scanpy as sc

# Load gene expression data
adata = sc.read_h5ad("adata.h5ad")
print(adata)

# With TCR/BCR data:
import muon as mu
mdata = mu.read("mudata.h5mu")  # Run load_data.py first to create this
```

## Quick Start (R)

```r
library(Seurat)
obj <- readRDS("../FINAL_integrated_object.rds")
```

## Notes

%s

Generated: %s
',
                    if (!is.null(airr_path)) "- `airr_rearrangement.tsv` - TCR/BCR in AIRR format\n- `mudata.h5mu` - Combined GEX+AIRR (created by load_data.py)" else "",
                    if (chunked) "- This dataset was exported using chunked HDF5 writing due to its large size\n- The matrix has >2 billion non-zero entries" else "",
                    Sys.time()
  )
  
  writeLines(readme, file.path(export_dir, "README.md"))
}


# =============================================================================
# AIRR Parsing Functions
# =============================================================================

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


parse_gene_string_to_airr <- function(cell_id, gene_str, cdr3_aa = NA, cdr3_nt = NA) {
  genes <- strsplit(gene_str, "\\.")[[1]]
  
  if (length(genes) == 0 || all(genes == "None") || all(genes == "NA")) return(NULL)
  
  airr_row <- list(
    cell_id = cell_id,
    sequence_id = paste0(cell_id, "_", gene_str),
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