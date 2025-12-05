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
  
  # --- 2. Export Metadata (without VDJ - that goes to AIRR) ---
  log_message("Exporting Metadata...")
  obs <- obj@meta.data
  obs$barcode <- rownames(obs)
  
  # Keep VDJ columns in obs for reference, but mark them
  vdj_cols <- c("CTaa", "CTnt", "CTgene", "CTstrict", 
                "clonalFrequency", "clonalProportion", "cloneSize",
                "TCR1", "TCR2", "cdr3_aa1", "cdr3_aa2", "cdr3_nt1", "cdr3_nt2")
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
  
  # --- 5. Export VDJ in AIRR format for scirpy ---
  log_message("Exporting VDJ data in AIRR format...")
  
  # Check if we have clonotype data attached
  has_vdj <- any(c("CTaa", "CTgene") %in% colnames(obj@meta.data))
  
  if (has_vdj) {
    log_message("Parsing VDJ from metadata into AIRR format...")
    airr_df <- parse_screpertoire_to_airr(obj@meta.data)
    
    if (!is.null(airr_df) && nrow(airr_df) > 0) {
      fwrite(airr_df, file = file.path(export_dir, "airr_rearrangement.tsv"), 
             sep = "\t", row.names = FALSE)
      log_message(paste0("Exported ", nrow(airr_df), " AIRR chains"))
    } else {
      log_message("Warning: Could not parse VDJ data into AIRR format")
    }
  }
  
  # --- 6. Generate Python Assembler Script ---
  log_message("Generating Python assembly script...")
  py_script_path <- file.path(export_dir, "assemble_mudata.py")
  
  py_code <- '#!/usr/bin/env python3
"""
Assembles exported R data into MuData format compatible with scirpy.
Requires: scanpy, anndata, muon, scirpy, pandas, scipy
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import sys
import warnings

# Check imports
try:
    import anndata
    import muon as mu
    import scirpy as ir
except ImportError as e:
    print(f"Missing required package: {e}")
    print("Install with: pip install muon scirpy")
    sys.exit(1)

print(f"Python Version: {sys.version}")
print(f"Pandas Version: {pd.__version__}")
print(f"Anndata Version: {anndata.__version__}")
print(f"Scirpy Version: {ir.__version__}")

print("\\n--- Assembling MuData from R export ---")

# 1. Load Gene Expression Data
print("Loading sparse matrix...")
try:
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

# Create GEX AnnData
adata_gex = anndata.AnnData(X=X, obs=obs.copy(), var=var)

# 2. Load Embeddings
print("Loading embeddings...")
for f in os.listdir("."):
    if f.startswith("obsm_") and f.endswith(".csv"):
        key = f.replace("obsm_", "").replace(".csv", "")
        emb = pd.read_csv(f, index_col=0).values
        adata_gex.obsm["X_" + key] = emb

# 3. Load AIRR data for scirpy
print("\\nProcessing VDJ/AIRR data for Scirpy...")

adata_airr = None

# Option A: Load from AIRR TSV if exported
if os.path.exists("airr_rearrangement.tsv"):
    print("Found AIRR rearrangement file, loading with scirpy...")
    try:
        adata_airr = ir.io.read_airr("airr_rearrangement.tsv")
        print(f"Loaded {adata_airr.n_obs} cells with AIRR data")
    except Exception as e:
        print(f"Warning: Could not load AIRR file: {e}")
        adata_airr = None

# Option B: Parse from scRepertoire columns in obs if AIRR file not available
if adata_airr is None and any(col in obs.columns for col in ["CTgene", "TCR1", "IGH"]):
    print("Parsing VDJ from scRepertoire columns...")
    try:
        adata_airr = parse_screpertoire_obs_to_airr(obs)
    except Exception as e:
        print(f"Warning: Could not parse scRepertoire data: {e}")

# 4. Create MuData object
print("\\nCreating MuData object...")

if adata_airr is not None and adata_airr.n_obs > 0:
    # Align indices - only keep cells present in both
    common_cells = adata_gex.obs_names.intersection(adata_airr.obs_names)
    print(f"Cells with GEX: {adata_gex.n_obs}")
    print(f"Cells with AIRR: {adata_airr.n_obs}")
    print(f"Cells in common: {len(common_cells)}")
    
    # Create MuData with both modalities
    mdata = mu.MuData({"gex": adata_gex, "airr": adata_airr})
    
    # Run scirpy preprocessing
    print("\\nRunning scirpy preprocessing...")
    try:
        ir.pp.index_chains(mdata)
        ir.tl.chain_qc(mdata)
        print("Chain QC complete!")
        
        # Show receptor type distribution
        if "receptor_type" in mdata["airr"].obs.columns:
            print("\\nReceptor type distribution:")
            print(mdata["airr"].obs["receptor_type"].value_counts())
    except Exception as e:
        print(f"Warning: scirpy preprocessing failed: {e}")
    
    # Save as MuData
    out_file = "final_export.h5mu"
    mdata.write(out_file)
    print(f"\\nSuccess! MuData saved to {out_file}")
    
    # Also save GEX-only h5ad for basic scanpy use
    adata_gex.write("gex_only.h5ad")
    print("GEX-only file saved to gex_only.h5ad")
    
else:
    # No AIRR data - save GEX only
    print("No AIRR data found, saving GEX-only AnnData...")
    out_file = "final_seurat_export.h5ad"
    adata_gex.write(out_file)
    print(f"Success! File saved to {out_file}")


def parse_screpertoire_obs_to_airr(obs_df):
    """
    Parse scRepertoire columns from obs into scirpy-compatible AIRR format.
    
    scRepertoire stores chains as:
    - CTgene: "TRAV12-1.TRAJ9.TRAC_TRBV5-1.TRBJ2-7.TRBC2" (underscore-separated)
    - CTaa: "CAVRGEGFQKLVF_CASSLTDRTYEQYF" (underscore-separated CDR3 aa)
    - CTnt: nucleotide sequences (underscore-separated)
    - TCR1/TCR2 or IGH/IGLC: individual chain gene info
    - cdr3_aa1/cdr3_aa2, cdr3_nt1/cdr3_nt2: individual CDR3 sequences
    """
    airr_cells = []
    
    for cell_id, row in obs_df.iterrows():
        cell = ir.io.AirrCell(cell_id=str(cell_id))
        
        # Determine if TCR or BCR based on available columns
        is_tcr = any(col in obs_df.columns for col in ["TCR1", "TCR2"])
        is_bcr = any(col in obs_df.columns for col in ["IGH", "IGLC", "IGK", "IGL"])
        
        chains_added = False
        
        if is_tcr:
            # Parse TCR chains
            for chain_idx, (chain_col, cdr3_aa_col, cdr3_nt_col) in enumerate([
                ("TCR1", "cdr3_aa1", "cdr3_nt1"),
                ("TCR2", "cdr3_aa2", "cdr3_nt2")
            ]):
                if chain_col in obs_df.columns and pd.notna(row.get(chain_col)):
                    chain_dict = ir.io.AirrCell.empty_chain_dict()
                    gene_str = str(row[chain_col])
                    
                    # Parse gene string like "TRAV12-1.TRAJ9.TRAC"
                    genes = gene_str.split(".")
                    
                    # Determine locus from V gene
                    if genes and len(genes) > 0:
                        v_gene = genes[0] if genes[0] != "None" else None
                        if v_gene:
                            if v_gene.startswith("TRA"):
                                chain_dict["locus"] = "TRA"
                            elif v_gene.startswith("TRB"):
                                chain_dict["locus"] = "TRB"
                            elif v_gene.startswith("TRG"):
                                chain_dict["locus"] = "TRG"
                            elif v_gene.startswith("TRD"):
                                chain_dict["locus"] = "TRD"
                            
                            chain_dict["v_call"] = v_gene
                        
                        # Parse other genes
                        for gene in genes[1:]:
                            if gene and gene != "None":
                                if "J" in gene[:3]:
                                    chain_dict["j_call"] = gene
                                elif "D" in gene[:3]:
                                    chain_dict["d_call"] = gene
                                elif "C" in gene[:3]:
                                    chain_dict["c_call"] = gene
                    
                    # Add CDR3 sequences
                    if cdr3_aa_col in obs_df.columns and pd.notna(row.get(cdr3_aa_col)):
                        chain_dict["junction_aa"] = str(row[cdr3_aa_col])
                    if cdr3_nt_col in obs_df.columns and pd.notna(row.get(cdr3_nt_col)):
                        chain_dict["junction"] = str(row[cdr3_nt_col])
                    
                    chain_dict["productive"] = True  # Assume productive if in scRepertoire output
                    
                    if chain_dict.get("locus"):
                        cell.add_chain(chain_dict)
                        chains_added = True
        
        elif is_bcr:
            # Parse BCR chains (IGH, IGK/IGL)
            for chain_col, cdr3_aa_col, cdr3_nt_col, locus in [
                ("IGH", "cdr3_aa1", "cdr3_nt1", "IGH"),
                ("IGLC", "cdr3_aa2", "cdr3_nt2", None),  # Could be IGK or IGL
            ]:
                if chain_col in obs_df.columns and pd.notna(row.get(chain_col)):
                    chain_dict = ir.io.AirrCell.empty_chain_dict()
                    gene_str = str(row[chain_col])
                    genes = gene_str.split(".")
                    
                    if genes and len(genes) > 0:
                        v_gene = genes[0] if genes[0] != "None" else None
                        if v_gene:
                            # Determine locus from V gene for light chains
                            if v_gene.startswith("IGH"):
                                chain_dict["locus"] = "IGH"
                            elif v_gene.startswith("IGK"):
                                chain_dict["locus"] = "IGK"
                            elif v_gene.startswith("IGL"):
                                chain_dict["locus"] = "IGL"
                            elif locus:
                                chain_dict["locus"] = locus
                            
                            chain_dict["v_call"] = v_gene
                        
                        for gene in genes[1:]:
                            if gene and gene != "None":
                                if "J" in gene[:4]:
                                    chain_dict["j_call"] = gene
                                elif "D" in gene[:4]:
                                    chain_dict["d_call"] = gene
                                elif "C" in gene[:4]:
                                    chain_dict["c_call"] = gene
                    
                    if cdr3_aa_col in obs_df.columns and pd.notna(row.get(cdr3_aa_col)):
                        chain_dict["junction_aa"] = str(row[cdr3_aa_col])
                    if cdr3_nt_col in obs_df.columns and pd.notna(row.get(cdr3_nt_col)):
                        chain_dict["junction"] = str(row[cdr3_nt_col])
                    
                    chain_dict["productive"] = True
                    
                    if chain_dict.get("locus"):
                        cell.add_chain(chain_dict)
                        chains_added = True
        
        # Fallback: parse from combined CT columns
        if not chains_added and "CTgene" in obs_df.columns and pd.notna(row.get("CTgene")):
            ct_gene = str(row["CTgene"])
            ct_aa = str(row.get("CTaa", "")) if pd.notna(row.get("CTaa")) else ""
            ct_nt = str(row.get("CTnt", "")) if pd.notna(row.get("CTnt")) else ""
            
            # Split by underscore (scRepertoire convention)
            gene_parts = ct_gene.split("_")
            aa_parts = ct_aa.split("_") if ct_aa else [""] * len(gene_parts)
            nt_parts = ct_nt.split("_") if ct_nt else [""] * len(gene_parts)
            
            for i, gene_str in enumerate(gene_parts):
                if gene_str and gene_str != "NA":
                    chain_dict = ir.io.AirrCell.empty_chain_dict()
                    genes = gene_str.split(".")
                    
                    if genes:
                        v_gene = genes[0] if genes[0] != "None" else None
                        if v_gene:
                            # Infer locus
                            if v_gene.startswith("TRA"):
                                chain_dict["locus"] = "TRA"
                            elif v_gene.startswith("TRB"):
                                chain_dict["locus"] = "TRB"
                            elif v_gene.startswith("TRG"):
                                chain_dict["locus"] = "TRG"
                            elif v_gene.startswith("TRD"):
                                chain_dict["locus"] = "TRD"
                            elif v_gene.startswith("IGH"):
                                chain_dict["locus"] = "IGH"
                            elif v_gene.startswith("IGK"):
                                chain_dict["locus"] = "IGK"
                            elif v_gene.startswith("IGL"):
                                chain_dict["locus"] = "IGL"
                            
                            chain_dict["v_call"] = v_gene
                        
                        for gene in genes[1:]:
                            if gene and gene != "None":
                                if "J" in gene[:4]:
                                    chain_dict["j_call"] = gene
                                elif "D" in gene[:4]:
                                    chain_dict["d_call"] = gene
                                elif "C" in gene[:4]:
                                    chain_dict["c_call"] = gene
                    
                    if i < len(aa_parts) and aa_parts[i] and aa_parts[i] != "NA":
                        chain_dict["junction_aa"] = aa_parts[i]
                    if i < len(nt_parts) and nt_parts[i] and nt_parts[i] != "NA":
                        chain_dict["junction"] = nt_parts[i]
                    
                    chain_dict["productive"] = True
                    
                    if chain_dict.get("locus"):
                        cell.add_chain(chain_dict)
                        chains_added = True
        
        if chains_added:
            airr_cells.append(cell)
    
    if airr_cells:
        return ir.io.from_airr_cells(airr_cells)
    return None


if __name__ == "__main__":
    # Change to export directory if running standalone
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=".", help="Export directory")
    args = parser.parse_args()
    
    os.chdir(args.dir)
    
    # Run assembly
    exec(open("assemble_mudata.py").read().split("if __name__")[0])
'
  
  writeLines(py_code, py_script_path)
  
  log_message("Export complete!")
  log_message(paste0("Run: cd ", export_dir, " && python3 assemble_mudata.py"))
}


#' Parse scRepertoire metadata columns into AIRR-compliant format
#' 
#' @param meta_data Seurat metadata data.frame with scRepertoire columns
#' @return data.frame in AIRR rearrangement format
parse_screpertoire_to_airr <- function(meta_data) {
  
  # Check what columns we have (Global check)
  has_individual_chains <- any(c("TCR1", "TCR2", "IGH", "IGLC") %in% colnames(meta_data))
  has_combined <- "CTgene" %in% colnames(meta_data)
  
  if (!has_individual_chains && !has_combined) {
    return(NULL)
  }
  
  # Use lapply to process each row independently
  # This returns a nested list: List of Cells -> List of Chains
  nested_rows <- lapply(seq_len(nrow(meta_data)), function(i) {
    cell_id <- rownames(meta_data)[i]
    row <- meta_data[i, , drop = FALSE]
    
    # Store rows for THIS cell only
    cell_airr_rows <- list()
    
    # 1. Try individual chain columns
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
            if (!is.null(airr_row)) {
              cell_airr_rows[[length(cell_airr_rows) + 1]] <- airr_row
            }
          }
        }
      }
    }
    
    # 2. Fallback to combined columns
    # We check length(cell_airr_rows) which is local to this specific cell iteration
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
            if (!is.null(airr_row)) {
              cell_airr_rows[[length(cell_airr_rows) + 1]] <- airr_row
            }
          }
        }
      }
    }
    
    return(cell_airr_rows)
  })
  
  # Flatten the list of lists (Cell -> Chains) into a single list of Chains
  # recursive = FALSE is critical here to only unwrap the top level (cells)
  airr_rows <- unlist(nested_rows, recursive = FALSE)
  
  if (!is.null(airr_rows) && length(airr_rows) > 0) {
    return(do.call(rbind, lapply(airr_rows, as.data.frame)))
  }
  return(NULL)
}
