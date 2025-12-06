# Seurat ↔ Scanpy/Scirpy Interoperability Guide

This guide enables users to seamlessly work with the same single-cell data in either R (Seurat/scRepertoire) or Python (Scanpy/Scirpy), achieving identical results on the same underlying data.

## Table of Contents

1. [Quick Start](#quick-start)
2. [File Formats](#file-formats)
3. [Loading Data](#loading-data)
4. [Field Equivalence](#field-equivalence)
5. [Troubleshooting](#troubleshooting)

---

## Quick Start

### For R Users (Starting Point)

The pipeline produces a Seurat object that can be used directly:

```r
library(Seurat)
library(scRepertoire)

# Load the integrated object
obj <- readRDS("results/FINAL_integrated_object.rds")

# Basic analysis
DimPlot(obj, group.by = "leiden.cluster")
FeaturePlot(obj, features = c("CD3E", "CD4", "CD8A"))

# VDJ analysis (if TCR/BCR data present)
clonalQuant(obj, cloneCall = "aa")
```

### For Python Users

The export creates scanpy/scirpy-compatible files:

```python
import scanpy as sc
import scirpy as ir
import muon as mu

# Option 1: Gene expression only
adata = sc.read_h5ad("results/09_scanpy_export/adata.h5ad")

# Option 2: Gene expression + TCR/BCR (recommended)
mdata = mu.read("results/09_scanpy_export/mudata.h5mu")

# Basic analysis
sc.pl.umap(mdata["gex"], color="leiden.cluster")

# VDJ analysis
ir.tl.define_clonotypes(mdata)
ir.pl.clonal_expansion(mdata)
```

---

## File Formats

| File | Format | Contains | Load With |
|------|--------|----------|-----------|
| `FINAL_integrated_object.rds` | R/Seurat | Everything | `readRDS()` |
| `adata.h5ad` | H5AD | GEX + metadata + embeddings | `sc.read_h5ad()` |
| `mudata.h5mu` | H5MU | GEX + AIRR (scirpy-ready) | `mu.read()` |
| `airr_rearrangement.tsv` | AIRR TSV | TCR/BCR chains | `ir.io.read_airr()` |

### When to Use Each Format

- **Seurat RDS**: Full R analysis, BPCells-backed for large data
- **H5AD**: Scanpy analysis, no VDJ needed
- **H5MU (MuData)**: Combined GEX + VDJ analysis with scirpy
- **AIRR TSV**: Standalone VDJ analysis or custom pipelines

---

## Loading Data

### R (Seurat)

```r
library(Seurat)

# Load object
obj <- readRDS("results/FINAL_integrated_object.rds")

# Check structure
print(obj)
# An object of class Seurat 
# 35000 features across 150000 samples within 1 assay 
# Active assay: RNA (35000 features, 3000 variable features)
# 3 dimensional reductions: pca, harmony, umap

# Access components
counts <- obj@assays$RNA@counts          # Count matrix
meta <- obj@meta.data                     # Cell metadata
umap <- Embeddings(obj, "umap")           # UMAP coordinates
```

### Python (Scanpy)

```python
import scanpy as sc

# Load H5AD
adata = sc.read_h5ad("results/09_scanpy_export/adata.h5ad")

# Check structure
print(adata)
# AnnData object with n_obs × n_vars = 150000 × 35000
#     obs: 'sample', 'cluster', 'nCount_RNA', ...
#     var: 'gene_symbols'
#     obsm: 'X_pca', 'X_harmony', 'X_umap'

# Access components
counts = adata.X                          # Count matrix (sparse)
meta = adata.obs                          # Cell metadata
umap = adata.obsm["X_umap"]               # UMAP coordinates
```

### Python (MuData with VDJ)

```python
import muon as mu

# Load MuData
mdata = mu.read("results/09_scanpy_export/mudata.h5mu")

# Check structure
print(mdata)
# MuData object with n_obs × n_vars = 150000 × 35000
#   2 modalities
#     gex: 150000 x 35000
#     airr: 90000 x 0

# Access modalities
gex = mdata["gex"]     # Gene expression (AnnData)
airr = mdata["airr"]   # Receptor data (AnnData)
```

---

## Field Equivalence

### Core Data Structures

| Concept | R (Seurat) | Python (Scanpy) |
|---------|------------|-----------------|
| Object | `Seurat` | `AnnData` |
| Cell count | `ncol(obj)` | `adata.n_obs` |
| Gene count | `nrow(obj)` | `adata.n_vars` |
| Cell names | `colnames(obj)` | `adata.obs_names` |
| Gene names | `rownames(obj)` | `adata.var_names` |

### Matrix Access

| Data | R (Seurat) | Python (Scanpy) |
|------|------------|-----------------|
| Raw counts | `obj@assays$RNA@counts` | `adata.X` or `adata.raw.X` |
| Normalized | `obj@assays$RNA@data` | `adata.layers["normalized"]` |
| Scaled | `obj@assays$RNA@scale.data` | `adata.layers["scaled"]` |

### Metadata

| Data | R (Seurat) | Python (Scanpy) |
|------|------------|-----------------|
| All metadata | `obj@meta.data` | `adata.obs` |
| Single column | `obj$cluster` | `adata.obs["cluster"]` |
| Gene metadata | `obj@assays$RNA@meta.features` | `adata.var` |

### Embeddings

| Reduction | R (Seurat) | Python (Scanpy) |
|-----------|------------|-----------------|
| PCA | `Embeddings(obj, "pca")` | `adata.obsm["X_pca"]` |
| UMAP | `Embeddings(obj, "umap")` | `adata.obsm["X_umap"]` |
| Harmony | `Embeddings(obj, "harmony")` | `adata.obsm["X_harmony"]` |

### VDJ/AIRR Data

| Data | R (scRepertoire) | Python (scirpy) |
|------|------------------|-----------------|
| TCR alpha V gene | `obj$TCR1` | `ir.get.airr(mdata, "v_call", "VJ_1")` |
| TCR beta V gene | `obj$TCR2` | `ir.get.airr(mdata, "v_call", "VDJ_1")` |
| CDR3 amino acid | `obj$CTaa` | `ir.get.airr(mdata, "junction_aa")` |
| Clonotype ID | `obj$CTstrict` | `mdata.obs["clone_id"]` |

---

## Troubleshooting

### H5AD Won't Load

```python
# Error: "Unable to open file"
# Solution: Check anndata version
pip install --upgrade anndata>=0.8.0

# Error: "KeyError: 'X'"  
# Solution: The matrix might be in a different location
import h5py
with h5py.File("adata.h5ad", "r") as f:
    print(list(f.keys()))  # Check structure
```

### Scirpy Chain QC Fails

```python
# Error: "No chain indices found"
# Solution: Run index_chains first
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)

# Error: "AttributeError: 'NoneType' object..."
# Solution: Ensure AIRR data was loaded correctly
print(mdata["airr"].obsm.keys())  # Should contain 'airr'
```

### Memory Issues

```python
# For very large datasets, use backed mode
adata = sc.read_h5ad("adata.h5ad", backed="r")

# Or load specific columns only
import pandas as pd
obs = pd.read_hdf("adata.h5ad", key="obs")
```

### Embedding Names Don't Match

```python
# R uses: pca, umap, harmony
# Python uses: X_pca, X_umap, X_harmony

# The export adds "X_" prefix automatically
# If you need to match R names:
adata.obsm["pca"] = adata.obsm["X_pca"]
```

### Missing VDJ Data in Python

```python
# If mudata.h5mu doesn't exist:
# 1. Check if AIRR file exists
import os
print(os.path.exists("airr_rearrangement.tsv"))

# 2. Create MuData manually
adata_gex = sc.read_h5ad("adata.h5ad")
adata_airr = ir.io.read_airr("airr_rearrangement.tsv")
mdata = mu.MuData({"gex": adata_gex, "airr": adata_airr})
ir.pp.index_chains(mdata)
mdata.write("mudata.h5mu")
```

---

## Version Compatibility

Tested with:

| Package | Version |
|---------|---------|
| R | ≥ 4.2.0 |
| Seurat | ≥ 5.0.0 |
| scRepertoire | ≥ 2.6.0 |
| Python | ≥ 3.9 |
| scanpy | ≥ 1.9.0 |
| anndata | ≥ 0.8.0 |
| scirpy | ≥ 0.13.0 |
| muon | ≥ 0.1.5 |

---

## Citation

If you use this interoperability pipeline, please cite:

- **Seurat**: Hao et al., Cell 2021
- **scRepertoire**: Borcherding et al., F1000Research 2020
- **Scanpy**: Wolf et al., Genome Biology 2018
- **Scirpy**: Sturm et al., Bioinformatics 2020
