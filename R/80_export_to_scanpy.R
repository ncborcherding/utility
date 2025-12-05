# R/80_export_to_scanpy.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(yaml)
  library(BPCells)
})

source("R/00_utils.R")

export_to_scanpy <- function(config_path = "config.yaml") {
  log_message("=== Starting Export to Scanpy/Scirpy ===")
  
  config <- yaml::read_yaml(config_path)
  input_path <- file.path(config$paths$results_dir, "FINAL_integrated_object.rds")
  
  export_dir <- file.path(config$paths$results_dir, "09_scanpy_export")
  
  if (!file.exists(input_path)) stop("Annotated object not found: ", input_path)
  obj <- readRDS(input_path)
  
  # --- 1. Export Counts (Sparse Matrix) ---
  log_message("Exporting RNA Counts...")
  counts <- obj@assays$RNA@layers$counts
  if (is.null(counts)) counts <- obj@assays$RNA@counts
  
  write_matrix_dir(
    mat = counts,
    dir = export_dir
  )
  
  # --- 2. Export Metadata ---
  log_message("Exporting Metadata...")
  obs <- obj@meta.data
  obs$barcode <- rownames(obs)
  
  # Ensure VDJ columns are character for safer CSV export
  vdj_cols <- c("CTaa", "CTnt", "CTgene", "CTstrict")
  for(c in vdj_cols) {
    if(c %in% colnames(obs)) obs[[c]] <- as.character(obs[[c]])
  }
  
  fwrite(obs, file = file.path(export_dir, "obs.csv"), row.names = FALSE)
  
  # --- 3. Export Features ---
  log_message("Exporting Features...")
  var <- data.frame(gene_symbols = rownames(obj), row.names = rownames(obj))
  fwrite(var, file = file.path(export_dir, "var.csv"), row.names = TRUE)
  
  # --- 4. Export Embeddings ---
  log_message("Exporting Embeddings...")
  for (red in names(obj@reductions)) {
    emb <- Embeddings(obj, reduction = red)
    write.csv(emb, file = file.path(export_dir, paste0("obsm_", red, ".csv")))
  }
  
  # --- 5. Generate Python Assembler Script ---
  # This script is specifically tuned for python=3.10, pandas>=2.0, and scirpy
  
  log_message("Generating Python assembly script...")
  py_script_path <- file.path(export_dir, "assemble_anndata.py")
  
  py_code <- sprintf('
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import anndata
import scirpy as ir
import sys

print(f"Python Version: {sys.version}")
print(f"Pandas Version: {pd.__version__}")
print(f"Anndata Version: {anndata.__version__}")

print("--- Assembling AnnData from R export ---")

# 1. Load Counts
print("Loading sparse matrix...")
try:
    # Scipy >= 1.10 mmread returns a coo_array or coo_matrix depending on version
    # We force CSR for AnnData efficiency
    X = sp.io.mmread("matrix.mtx.gz").tocsr().T
except Exception as e:
    print(f"Error reading matrix: {e}")
    sys.exit(1)

obs = pd.read_csv("obs.csv", index_col="barcode")
var = pd.read_csv("var.csv", index_col=0)

# Verify dimensions
if X.shape[0] != obs.shape[0]:
    print(f"Transposing matrix. X: {X.shape}, Obs: {obs.shape}")
    X = X.T

# 2. Create AnnData
adata = anndata.AnnData(X=X, obs=obs, var=var)

# 3. Load Embeddings
print("Loading embeddings...")
files = os.listdir(".")
for f in files:
    if f.startswith("obsm_") and f.endswith(".csv"):
        key = f.replace("obsm_", "").replace(".csv", "")
        # pandas 2.0 read_csv is standard
        emb = pd.read_csv(f, index_col=0).values
        adata.obsm["X_" + key] = emb

# 4. Processing VDJ for Scirpy (Reconstructing AIRR)
print("Processing VDJ data for Scirpy...")

if "CTaa" in adata.obs.columns and "CTgene" in adata.obs.columns:
    print("Found scRepertoire columns. Attempting to parse into IR format...")
    
    # scRepertoire stores chains combined (e.g., TRA:ACGT;TRB:ACGT)
    # We need to convert this to an object scirpy understands.
    # Strategy: Create a dummy AIRR table or use scirpy from_anndata logic if data is formatted.
    # Best approach for simple metadata: Let scirpy define clonotypes based on the columns we have.
    
    # However, to get FULL scirpy functionality (chain pairing plots), scirpy needs "ir" data.
    # We will rely on the existing clonotype ID from R if available, 
    # OR instruct scirpy to define clonotypes based on CTnt.
    
    try:
        # map R clonotype info to scirpy if needed
        if "CTstrict" in adata.obs.columns:
             ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only", key_added="clonotype")
    except Exception as e:
        print(f"Auto-clonotype definition skipped: {e}")
        print("Note: VDJ metadata is preserved in .obs for manual scirpy ingestion.")

# 5. Write h5ad
out_file = "final_seurat_export.h5ad"
adata.write(out_file)
print(f"Success! File saved to {out_file}")
', export_dir)
  
  writeLines(py_code, py_script_path)
  log_message("Export complete. Run 'python3 assemble_anndata.py' in ", export_dir)
}