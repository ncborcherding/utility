#!/usr/bin/env python3
"""Verify H5AD export with int64 index check."""
import sys
try:
    import scanpy as sc
    import numpy as np
except ImportError as e:
    print(f"Missing: {e}")
    sys.exit(1)

print(f"scanpy: {sc.__version__}")
print("\nLoading adata.h5ad...")

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
    
    print("\n✓ Verification passed!")
except Exception as e:
    print(f"✗ {e}")
    sys.exit(1)

