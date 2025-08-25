# py/scanvi_integration.py
#
# Python script to run scANVI integration, called from R via reticulate.

import scvi
import scanpy as sc
import anndata as ad
import yaml
import logging
import pandas as pd

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(message)s')

def run_scanvi_integration(h5ad_path: str, config_path: str):
    """
    Sets up and trains a scANVI model on the given AnnData object.

    Args:
        h5ad_path (str): Path to the .h5ad file containing the single-cell data.
        config_path (str): Path to the YAML configuration file.
    """

    # --- 1. Read Config and AnnData ---
    logging.info(f"Reading configuration from: {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    scanvi_config = config['methods']['scanvi']
    batch_key = config['methods']['harmony']['group_by_vars'][0]
    labels_key = scanvi_config['labels_key']
    unlabeled_category = scanvi_config['unlabeled_category']

    logging.info(f"Reading AnnData object from: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)

    # Ensure the labels column is of a compatible type (e.g., string/category)
    if labels_key in adata.obs:
        if pd.api.types.is_numeric_dtype(adata.obs[labels_key]):
             adata.obs[labels_key] = adata.obs[labels_key].astype(str)
        adata.obs[labels_key] = adata.obs[labels_key].astype("category")
    else:
        raise ValueError(f"Labels key '{labels_key}' not found in adata.obs.")


    # --- 2. Setup and Train scANVI Model ---
    logging.info("Setting up scANVI model...")
    scvi.model.SCANVI.setup_anndata(
        adata,
        batch_key=batch_key,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category
    )

    # Train from scratch (not from a pre-trained scVI model) for simplicity
    model = scvi.model.SCANVI(
        adata,
        n_layers=scanvi_config['n_layers'],
        n_latent=scanvi_config['n_latent'],
        n_hidden=scanvi_config['n_hidden']
    )

    logging.info(f"Training scANVI model for {scanvi_config['max_epochs']} epochs...")
    # scANVI training is often shorter than scVI
    model.train(
        max_epochs=scanvi_config['max_epochs'],
        early_stopping=True,
        lr=scanvi_config['learning_rate'],
        batch_size=scanvi_config['batch_size']
    )

    # --- 3. Get Latent Representation and Save ---
    logging.info("Training complete. Getting latent representation...")
    latent_embedding = model.get_latent_representation()

    # Store the embedding in the AnnData object
    adata.obsm["X_scanvi"] = latent_embedding

    logging.info(f"Saving updated AnnData object back to: {h5ad_path}")
    adata.write_h5ad(h5ad_path, compression="gzip")

    logging.info("scANVI integration script finished successfully.")

# This script is designed to be called from R.
