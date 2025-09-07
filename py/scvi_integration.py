# py/scvi_integration.py
#
# Python script to run scVI integration, called from R via reticulate.

import logging
import inspect
import yaml
import anndata as ad
import numpy as np
import scipy.sparse as sp
import pandas as pd  # <-- needed for categorical handling

import scvi

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(message)s')
log = logging.getLogger(__name__)


# ---------------------- Helpers ---------------------- #

def _is_integerish_sparse_or_array(X, atol=1e-6):
    """Return True if X values are (near) integers."""
    data = X.data if sp.issparse(X) else np.asarray(X)
    return np.allclose(data, np.round(data), atol=atol)


def _ensure_counts_in_X(adata):
    """
    Ensure raw integer count data are stored in:
      - adata.X (CSR, int32), and
      - adata.layers['counts'] (same as X)
    Preserve any normalized/float X as layers['X_norm'] for reference.
    """
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X

    if not _is_integerish_sparse_or_array(counts):
        log.warning("Counts appear non-integer; rounding to nearest integers for scVI.")
        if sp.issparse(counts):
            counts = counts.tocsr(copy=True)
            counts.data = np.rint(counts.data)
        else:
            counts = sp.csr_matrix(np.rint(np.asarray(counts)))

    # Preserve original X if it looks normalized / float
    if "X_norm" not in adata.layers and not _is_integerish_sparse_or_array(adata.X):
        adata.layers["X_norm"] = adata.X.copy()

    # Force X -> CSR int32
    counts_csr = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    counts_csr.data = counts_csr.data.astype(np.int32, copy=False)
    adata.X = counts_csr
    adata.layers["counts"] = adata.X  # mirror


def _collapse_small_batches(adata, batch_key, min_cells=3):
    """
    Map any batch categories with < min_cells to a single '<batch_key>_other' bucket.
    Works for both object dtype and pandas.Categorical dtype.
    """
    if batch_key is None or batch_key not in adata.obs:
        return

    ser = adata.obs[batch_key]
    other_label = f"{batch_key}_other"

    # Value counts on the *actual values* (works for categoricals too)
    vc = ser.value_counts(dropna=False)
    small = vc[vc < min_cells].index

    if len(small) == 0:
        return

    # If categorical, add the new category first, then replace; optionally prune old small categories
    if pd.api.types.is_categorical_dtype(ser):
        # add new category if not present
        if other_label not in ser.cat.categories:
            ser = ser.cat.add_categories([other_label])

        # replace small categories with 'other'
        ser = ser.where(~ser.isin(small), other_label)

        # (optional) remove the now-unused small categories if they still exist
        # Only remove categories that are actually present in the category set
        to_remove = [c for c in small if c in list(ser.cat.categories)]
        if to_remove:
            try:
                ser = ser.cat.remove_categories(to_remove)
            except Exception:
                # If removal fails for any reason, keep categories; not critical for training
                pass

        adata.obs[batch_key] = ser
    else:
        # Object/string dtype
        adata.obs[batch_key] = ser.where(~ser.isin(small), other_label)


def _pick_accelerator():
    """Prefer Apple Silicon 'mps' when available; otherwise cpu."""
    try:
        import torch
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
    except Exception:
        pass
    return "cpu"


def _maybe_kwarg(sig, name):
    """Return True if a callable signature supports a kwarg."""
    return name in sig.parameters


# ---------------------- Main entry ---------------------- #

def run_scvi_integration(h5ad_path: str, config_path: str):
    """
    Sets up and trains an scVI model on the given AnnData object.

    Args:
        h5ad_path (str): Path to the .h5ad file containing the single-cell data.
        config_path (str): Path to the YAML configuration file.
    """
    # --- 1) Read config and data ---
    log.info("Reading configuration from: %s", config_path)
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    scvi_cfg = cfg.get("methods", {}).get("scvi", {}) or {}
    # keep harmony's group_by_vars for consistency
    batch_key = (cfg.get("methods", {})
                   .get("harmony", {})
                   .get("group_by_vars", [None]))[0]

    log.info("Reading AnnData object from: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)

    # --- 2) Sanity & preprocessing for scVI ---
    _ensure_counts_in_X(adata)
    if batch_key is not None and batch_key in adata.obs:
        _collapse_small_batches(adata, batch_key, min_cells=3)
    else:
        batch_key = None

    # Ensure setup uses counts layer
    log.info("Setting up scVI model...")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)

    n_layers = int(scvi_cfg.get("n_layers", 2))
    n_latent = int(scvi_cfg.get("n_latent", 30))
    n_hidden = int(scvi_cfg.get("n_hidden", 128))
    model = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent, n_hidden=n_hidden)

    # --- 3) Train (version-robust kwargs handling) ---
    max_epochs = int(scvi_cfg.get("max_epochs", 400))
    lr_val = float(scvi_cfg.get("lr", scvi_cfg.get("learning_rate", 1e-3)))
    batch_size = scvi_cfg.get("batch_size", None)

    train_sig = inspect.signature(model.train)
    trainer_kwargs = dict(
        max_epochs=max_epochs,
        enable_checkpointing=False,
        check_val_every_n_epoch=5,
        early_stopping=True,
        accelerator=_pick_accelerator(),
        devices="auto",
    )

    # Include batch_size only if supported by this scvi-tools version
    if _maybe_kwarg(train_sig, "batch_size") and batch_size is not None:
        trainer_kwargs["batch_size"] = int(batch_size)

    # Route LR through TrainingPlan only if supported
    if _maybe_kwarg(train_sig, "plan_kwargs"):
        from scvi.train import TrainingPlan
        plan_sig = inspect.signature(TrainingPlan.__init__)
        if _maybe_kwarg(plan_sig, "lr"):
            trainer_kwargs["plan_kwargs"] = {"lr": lr_val}
        elif _maybe_kwarg(plan_sig, "learning_rate"):
            trainer_kwargs["plan_kwargs"] = {"learning_rate": lr_val}
        # else: no explicit LR override supported; fall back to defaults

    log.info("Training scVI model for %d epochs...", max_epochs)
    model.train(**trainer_kwargs)

    # --- 4) Save latent representation back into AnnData ---
    log.info("Writing X_scvi to obsm and saving AnnData...")
    adata.obsm["X_scvi"] = model.get_latent_representation()
    adata.write_h5ad(h5ad_path, compression="gzip")
    log.info("Done.")
