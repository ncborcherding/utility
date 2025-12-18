#!/usr/bin/env python3
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
    airr_df = pd.read_csv(airr_file, sep="\t", compression="gzip")
    temp_tsv = tempfile.mktemp(suffix=".tsv")
    airr_df.to_csv(temp_tsv, sep="\t", index=False)
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

