# py/scvi_integration.py
# -----------------------------------------------------------------------------
# Memory-lean scVI integration that can respect Seurat HVGs stamped in var['highly_variable'].
# -----------------------------------------------------------------------------

import os
import gc
import logging
import inspect
import yaml
import anndata as ad
import numpy as np
import scipy.sparse as sp
import pandas as pd

# --- Safety knobs for macOS / MPS and thread contention ---
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF",
    "expandable_segments:True,garbage_collection_threshold:0.9,max_split_size_mb:64")
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
os.environ.setdefault("PYTORCH_MPS_HIGH_WATERMARK_RATIO", "0.0")  # avoid aggressive caching on MPS
os.environ.setdefault("OMP_NUM_THREADS", "1")       # tame BLAS threads
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import torch  # import before scvi so we can set a couple of knobs
try:
    torch.set_num_threads(1)   # keep CPU threads in check when we run on CPU
except Exception:
    pass

import scvi   # now bring in scvi

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(message)s')
log = logging.getLogger(__name__)


# ---------------------- Helpers ---------------------- #

def _is_integerish_sparse_or_array(X, atol=1e-6):
    """Return True if X values are (near) integers."""
    data = X.data if sp.issparse(X) else np.asarray(X)
    return np.allclose(data, np.round(data), atol=atol)


def _ensure_counts_in_X(adata, prefer_layer="counts"):
    """
    Ensure raw integer count data are stored in adata.X (CSR) and adata.layers['counts'].
    - If source is float-ish, round to nearest integer.
    - Cast values to uint16 when possible, else int32 to reduce memory.
    - Preserve original X if it looks normalized in adata.layers['X_norm'].
    """
    src = adata.layers[prefer_layer] if prefer_layer in adata.layers else adata.X

    if not _is_integerish_sparse_or_array(src):
        log.warning("Counts appear non-integer; rounding to nearest integers for scVI.")
        if sp.issparse(src):
            src = src.tocsr(copy=True)
            src.data = np.rint(src.data)
        else:
            src = sp.csr_matrix(np.rint(np.asarray(src)))

    counts_csr = src.tocsr() if sp.issparse(src) else sp.csr_matrix(src)

    maxv = counts_csr.data.max() if counts_csr.nnz > 0 else 0
    if maxv <= np.iinfo(np.uint16).max:
        counts_csr.data = counts_csr.data.astype(np.uint16, copy=False)
    else:
        counts_csr.data = counts_csr.data.astype(np.int32, copy=False)

    # If the original X looked float-like (normalized), save it for provenance
    if "X_norm" not in adata.layers and not _is_integerish_sparse_or_array(adata.X):
        try:
            adata.layers["X_norm"] = adata.X.copy()
        except Exception:
            pass

    adata.X = counts_csr
    adata.layers["counts"] = adata.X  # mirror


def _collapse_small_batches(adata, batch_key, min_cells=3):
    """Collapse tiny batches to a single 'other' label to stabilize training."""
    if batch_key is None or batch_key not in adata.obs:
        return
    ser = adata.obs[batch_key]
    other_label = f"{batch_key}_other"
    vc = ser.value_counts(dropna=False)
    small = vc[vc < max(1, int(min_cells))].index
    if len(small) == 0:
        return

    if pd.api.types.is_categorical_dtype(ser):
        if other_label not in ser.cat.categories:
            ser = ser.cat.add_categories([other_label])
        ser = ser.where(~ser.isin(small), other_label)
        try:
            ser = ser.cat.remove_unused_categories()
        except Exception:
            pass
        adata.obs[batch_key] = ser
    else:
        adata.obs[batch_key] = ser.where(~ser.isin(small), other_label)


def _pick_accelerator():
    """auto-detect a good default accelerator."""
    acc = "cpu"
    try:
        import torch
        if torch.cuda.is_available():
            acc = "gpu"
        elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
            acc = "mps"
    except Exception:
        pass
    return acc


def _maybe_kwarg(sig, name):
    """Return True if a callable signature accepts the given keyword argument."""
    try:
        return name in sig.parameters
    except Exception:
        return False


def _prune_anndata(
    adata,
    keep_layers=("counts",),
    drop_raw=True,
    drop_obsm=True,
    drop_obsp=True,
    drop_varm=True,
    drop_layers_except=True,
):
    """Remove heavy, unused containers to reduce RAM."""
    if drop_layers_except:
        new_layers = {}
        for k in keep_layers:
            if k in adata.layers:
                new_layers[k] = adata.layers[k]
        adata.layers = new_layers

    if drop_raw and getattr(adata, "raw", None) is not None:
        adata.raw = None

    if drop_obsm and len(adata.obsm) > 0:
        adata.obsm.clear()
    if drop_obsp and len(adata.obsp) > 0:
        adata.obsp.clear()
    if drop_varm and len(adata.varm) > 0:
        adata.varm.clear()

    return adata


def _select_hvgs(adata, n_top=3000, flavor="seurat_v3", use_var_flag=True):
    """
    Subset to HVGs.
    Priority:
      1) If use_var_flag and var['highly_variable'] exists and has any True, subset by it.
      2) Else, if n_top > 0, compute HVGs (scanpy preferred).
      3) Else, do nothing.
    """
    if use_var_flag and "highly_variable" in adata.var.columns:
        hv = adata.var["highly_variable"].values
        if hv.dtype != bool:
            hv = hv.astype(bool)
        if hv.any():
            adata._inplace_subset_var(hv)
            log.info("Using precomputed HVGs from adata.var['highly_variable'] (n=%d).", int(hv.sum()))
            return adata
        else:
            log.info("var['highly_variable'] present but empty; will compute HVGs if requested.")

    if n_top is None or n_top <= 0:
        log.info("No HVG selection requested (keeping all genes).")
        return adata

    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if _HAS_SCANPY:
        log.info("Computing %d HVGs with scanpy (flavor=%s).", int(n_top), flavor)
        tmp = ad.AnnData(X=X, obs=adata.obs.copy(), var=adata.var.copy())
        sc.pp.highly_variable_genes(tmp, n_top_genes=int(n_top), flavor=str(flavor))
        hv = tmp.var["highly_variable"].values
        adata._inplace_subset_var(hv)
    else:
        log.info("scanpy not available; approximating HVGs via highest variance (n=%d).", int(n_top))
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        mean = np.asarray(X.mean(axis=0)).ravel()
        sq = X.copy()
        sq.data = sq.data.astype(np.float64) ** 2
        mean_sq = np.asarray(sq.mean(axis=0)).ravel()
        var = mean_sq - mean ** 2
        hv_idx = np.argsort(var)[::-1][:int(n_top)]
        hv = np.zeros(X.shape[1], dtype=bool)
        hv[hv_idx] = True
        adata._inplace_subset_var(hv)

    log.info("HVG selection complete; genes retained: %d.", adata.n_vars)
    return adata


# ---------------------- Main entry ---------------------- #

def run_scvi_integration(h5ad_path: str, config_path: str):
    """
    Train scVI on the given AnnData and write latent to obsm['X_scvi'] in-place.

    Args:
        h5ad_path (str): Path to input/output .h5ad file (will be updated).
        config_path (str): Path to YAML configuration with 'methods' and 'memory' sections.
    """
    import gc
    import inspect
    import anndata as ad
    import yaml
    import numpy as np

    log.info("Reading configuration from: %s", config_path)
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f) or {}

    methods  = cfg.get("methods", {}) or {}
    scvi_cfg = (methods.get("scvi", {}) or {}).copy()
    mem_cfg  = (cfg.get("memory", {}) or {}).copy()

    # Use the same batch variable as Harmony for consistency
    batch_key = (methods.get("harmony", {}).get("group_by_vars", [None]))[0]

    log.info("Reading AnnData object from: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)

    # ------------------ Preprocess & Memory Pruning ------------------ #
    _ensure_counts_in_X(adata)

    _collapse_small_batches(
        adata,
        batch_key if batch_key in adata.obs else None,
        min_cells=int(mem_cfg.get("min_cells_per_batch", 3)),
    )
    if batch_key not in adata.obs:
        batch_key = None

    adata = _prune_anndata(
        adata,
        keep_layers=tuple(mem_cfg.get("keep_layers", ["counts"])),
        drop_raw=bool(mem_cfg.get("drop_raw", True)),
        drop_obsm=bool(mem_cfg.get("drop_obsm", True)),
        drop_obsp=bool(mem_cfg.get("drop_obsp", True)),
        drop_varm=bool(mem_cfg.get("drop_varm", True)),
        drop_layers_except=bool(mem_cfg.get("drop_layers_except", True)),
    )

    # HVG policy:
    # - If you stamped Seurat HVGs into var['highly_variable'] and set
    #   use_var_highly_variable: true & hvg_n: null, we keep exactly that set.
    use_hvg_flag = bool(scvi_cfg.get("use_var_highly_variable", True))
    hvg_n_raw    = scvi_cfg.get("hvg_n", 3000)
    hvg_n        = 0 if hvg_n_raw is None else int(hvg_n_raw)
    adata = _select_hvgs(
        adata,
        n_top=hvg_n,
        flavor=str(scvi_cfg.get("hvg_flavor", "seurat_v3")),
        use_var_flag=use_hvg_flag,
    )

    gc.collect()

    # ------------------ scVI Setup ------------------ #
    import scvi

    log.info("Setting up scVI (layer='counts', batch_key=%s)...", batch_key)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)

    n_layers = int(scvi_cfg.get("n_layers", 2))
    n_latent = int(scvi_cfg.get("n_latent", 30))
    n_hidden = int(scvi_cfg.get("n_hidden", 128))

    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
    )

    # ------------------ Trainer Settings (safe defaults) ------------------ #
    max_epochs = int(scvi_cfg.get("max_epochs", 400))
    lr_val     = float(scvi_cfg.get("lr", scvi_cfg.get("learning_rate", 1e-3)))
    batch_size = scvi_cfg.get("batch_size", 512)

    # Dataloader kept conservative to avoid RAM spikes & R freezes
    dl_kwargs = {
        "num_workers": int(mem_cfg.get("num_workers", 0)),
        "pin_memory": bool(mem_cfg.get("pin_memory", False)),
        "persistent_workers": False,
    }

    train_sig = inspect.signature(model.train)
    trainer_kwargs = dict(
        max_epochs=max_epochs,
        enable_checkpointing=False,
        check_val_every_n_epoch=int(scvi_cfg.get("check_val_every_n_epoch", 5)),
        early_stopping=bool(scvi_cfg.get("early_stopping", True)),
        accelerator=scvi_cfg.get("accelerator", "auto"),
        devices=scvi_cfg.get("devices", "auto"),
        enable_progress_bar=bool(scvi_cfg.get("enable_progress_bar", True)),
        data_loader_kwargs=dl_kwargs if _maybe_kwarg(train_sig, "data_loader_kwargs") else None,
    )
    # Drop Nones
    trainer_kwargs = {k: v for k, v in trainer_kwargs.items() if v is not None}

    # Resolve 'auto' accelerator early
    if trainer_kwargs.get("accelerator") == "auto":
        trainer_kwargs["accelerator"] = _pick_accelerator()

    # Respect batch_size only if supported by this scvi-tools version
    if _maybe_kwarg(train_sig, "batch_size") and batch_size is not None:
        trainer_kwargs["batch_size"] = int(batch_size)

    # Precision handling (important for macOS/MPS)
    precision = scvi_cfg.get("precision", None)   # e.g., "16-mixed", "bf16-mixed", "32-true"
    if precision and _maybe_kwarg(train_sig, "precision"):
        trainer_kwargs["precision"] = precision

    # If MPS is selected, force full precision and modest batch sizes (AMP unsupported on MPS)
    if str(trainer_kwargs.get("accelerator")) == "mps":
        if _maybe_kwarg(train_sig, "precision"):
            if trainer_kwargs.get("precision") != "32-true":
                log.info("MPS detected: overriding precision -> 32-true")
                trainer_kwargs["precision"] = "32-true"
        if "batch_size" in trainer_kwargs and trainer_kwargs["batch_size"] > 128:
            log.info("MPS detected: capping batch_size -> 128")
            trainer_kwargs["batch_size"] = 128

    # If CPU was requested, ensure precision is sane
    if str(trainer_kwargs.get("accelerator")) == "cpu":
        if _maybe_kwarg(train_sig, "precision"):
            trainer_kwargs.setdefault("precision", "32-true")

    # Route learning rate through TrainingPlan when available
    if _maybe_kwarg(train_sig, "plan_kwargs"):
        try:
            from scvi.train import TrainingPlan
            plan_sig = inspect.signature(TrainingPlan.__init__)
            if _maybe_kwarg(plan_sig, "lr"):
                trainer_kwargs["plan_kwargs"] = {"lr": lr_val}
            elif _maybe_kwarg(plan_sig, "learning_rate"):
                trainer_kwargs["plan_kwargs"] = {"learning_rate": lr_val}
        except Exception:
            pass

    log.info(
        "Training scVI (epochs=%d, batch_size=%s, accelerator=%s, precision=%s)...",
        max_epochs,
        trainer_kwargs.get("batch_size", "auto"),
        trainer_kwargs.get("accelerator", "cpu"),
        trainer_kwargs.get("precision", "default"),
    )

    # ------------------ Train ------------------ #
    try:
        model.train(**trainer_kwargs)
    finally:
        # Encourage Python to free temporary tensors even if training fails
        gc.collect()

    # ------------------ Save Latent & Write Back ------------------ #
    log.info("Extracting latent representation to obsm['X_scvi'] and writing file...")
    adata.obsm["X_scvi"] = model.get_latent_representation()

    # Free model before write to shrink peak RSS
    del model
    gc.collect()

    adata.write_h5ad(h5ad_path, compression="gzip")
    log.info("Done writing to: %s", h5ad_path)
