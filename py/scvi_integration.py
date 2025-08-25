# py/scvi_integration.py
#
# Python script to run scVI integration, called from R via reticulate.

import scvi
import scanpy as sc
import anndata as ad
import yaml
import logging

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(message)s')

def run_scvi_integration(h5ad_path: str, config_path: str):
    """
    Sets up and trains an scVI model on the given AnnData object.

    Args:
        h5ad_path (str): Path to the .h5ad file containing the single-cell data.
        config_path (str): Path to the YAML configuration file.
    """

    # --- 1. Read Config and AnnData ---
    logging.info(f"Reading configuration from: {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    scvi_config = config['methods']['scvi']
    # Use the same batch variable as harmony for consistency
    batch_key = config['methods']['harmony']['group_by_vars'][0]

    logging.info(f"Reading AnnData object from: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)

    # --- 2. Setup and Train scVI Model ---
    logging.info("Setting up scVI model...")
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=batch_key
    )

    model = scvi.model.SCVI(
        adata,
        n_layers=scvi_config['n_layers'],
        n_latent=scvi_config['n_latent'],
        n_hidden=scvi_config['n_hidden']
    )

    logging.info(f"Training scVI model for {scvi_config['max_epochs']} epochs...")
    model.train(
        max_epochs=scvi_config['max_epochs'],
        early_stopping=True,
        lr=scvi_config['learning_rate'],
        batch_size=scvi_config['batch_size'],
        plan_kwargs={'lr_patience': 8, 'lr_factor': 0.1} # Add some robustness
    )

    # --- 3. Get Latent Representation and Save ---
    logging.info("Training complete. Getting latent representation...")
    latent_embedding = model.get_latent_representation()

    # Store the embedding in the AnnData object
    adata.obsm["X_scvi"] = latent_embedding

    logging.info(f"Saving updated AnnData object back to: {h5ad_path}")
    adata.write_h5ad(h5ad_path, compression="gzip")

    logging.info("scVI integration script finished successfully.")

# This script is designed to be called from R, so it doesn't need a __main__ block.
# The R script will source this file and call the run_scvi_integration function.
