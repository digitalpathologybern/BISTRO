#!/usr/bin/env python
"""
BISTRO HVG Selection Benchmarking
==================================

Processes a SINGLE (dataset, layer) pair. SLURM / Nextflow handles
parallelism across all combinations.

Benchmark phases (Section 2.4 / Section 4.8 of the manuscript):
  Phase 2a  Clustering stability   (bootstrap Leiden ARI)
  Phase 3   HVG vs. random control (ARI advantage over random gene sets)
  Phase 4a  Spatial coherence      (kNN label coherence in physical space, Eq. 8)
  Phase 4b  Marker gene AUROC      (Wilcoxon DE top-5 marker discriminability)

Usage
-----
    python hvg_benchmark.py \\
        --zarr /path/to/dataset.zarr \\
        --nextflow_output /path/to/nextflow_output/ \\
        --technology CosMx \\
        --layer cpm \\
        --output_dir /path/to/output/

    # List available layers
    python hvg_benchmark.py \\
        --zarr /path/to/dataset.zarr \\
        --nextflow_output /path/to/nextflow_output/ \\
        --technology CosMx \\
        --list_layers
"""

from __future__ import annotations
 
import argparse
import gc
import json
import os
import shutil
import sys
import tempfile
import time
import warnings
from itertools import combinations
from pathlib import Path
 
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.decomposition import PCA as skPCA
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
 
try:
    import igraph as ig
    import leidenalg
    HAS_LEIDENALG = True
except ImportError:
    HAS_LEIDENALG = False

warnings.filterwarnings("ignore")

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'utils'))
from helpers import parse_layer as _parse_layer
from data_loader import (
    load_spatialdata, discover_norm_layers, load_layer,
    get_dense_layer, get_spatial_coordinates,
)


# ============================================================================
# LOGGING
# Source: 03_HVG_benchmark.py, line 188
# ============================================================================

def log(msg: str) -> None:
    """Print a timestamped log line (immediately flushed)."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================================
# GENE SELECTION
# Source: 03_HVG_benchmark.py, lines 286-306
# ============================================================================

def load_hvg_mask(adata, layer_name: str, fraction: float,
                  nextflow_path: str) -> np.ndarray:
    """Load pre-computed boolean HVG mask from Nextflow CSV."""
    hvg_path = os.path.join(
        nextflow_path, "hvg",
        f"top_{fraction}_selected_hvg_per_method.csv",
    )
    hvg_df = pd.read_csv(hvg_path, index_col=0)
    selected_genes = hvg_df.index[hvg_df[layer_name] == True]  # noqa: E712
    return adata.var_names.isin(selected_genes)


def random_gene_mask(n_genes: int, fraction: float, rng) -> np.ndarray:
    """Return a random boolean mask selecting ``fraction`` of genes."""
    n_select = max(2, int(n_genes * fraction))
    if fraction >= 1.0:
        return np.ones(n_genes, dtype=bool)
    idx = rng.choice(n_genes, size=n_select, replace=False)
    mask = np.zeros(n_genes, dtype=bool)
    mask[idx] = True
    return mask


# ============================================================================
# FAST DOWNSTREAM: PCA + cKDTree/igraph Leiden
# ============================================================================

def _single_bootstrap(b, pca_embedding, n_cells, n_sub, cfg, base_seed):
    """Worker function for parallelized bootstrap iterations."""
    rng = np.random.default_rng(base_seed + b)
    cell_idx = np.sort(rng.choice(n_cells, size=n_sub, replace=False))

    try:
        labels = leiden_from_pca(
            pca_embedding[cell_idx],
            k=cfg["k_neighbors"],
            resolution=cfg["resolution"],
            seed=base_seed + b,
        )
        return b, cell_idx, labels
    except Exception:
        return b, None, None


def compute_pca_embedding(
    adata, layer_name: str, gene_mask: np.ndarray,
    n_components: int, seed: int = 42,
) -> np.ndarray:
    """
    Compute PCA embedding on HVG-filtered expression matrix.

    Handles NaN values from normalization methods by replacing with 0.

    Source: 03_HVG_benchmark.py, line 331.
    """
    X = get_dense_layer(adata[:, gene_mask], layer_name)

    nan_count = np.isnan(X).sum()
    if nan_count > 0:
        nan_pct = nan_count / X.size * 100
        log(f"      Warning: {nan_count} NaN values ({nan_pct:.2f}%) replaced with 0")
        np.nan_to_num(X, copy=False, nan=0.0)

    n_comps = min(n_components, X.shape[1] - 1, X.shape[0] - 1)
    pca = skPCA(n_components=n_comps, random_state=seed)
    embedding = pca.fit_transform(X).astype(np.float32)
    del X
    return embedding


def leiden_from_pca(
    pca_emb: np.ndarray, k: int, resolution: float, seed: int = 42,
) -> np.ndarray:
    """
    Build kNN graph from PCA embedding using cKDTree, construct an
    igraph object, and run Leiden. Returns integer cluster labels.
    """
    n_cells = pca_emb.shape[0]
    tree = cKDTree(pca_emb)
    distances, indices = tree.query(pca_emb, k=k + 1, workers=1)
    distances = distances[:, 1:]
    indices = indices[:, 1:]

    if HAS_LEIDENALG:
        sources = np.repeat(np.arange(n_cells), k)
        targets = indices.ravel()
        weights = np.exp(-distances.ravel())

        g = ig.Graph(n=n_cells, edges=list(zip(sources, targets)), directed=True)
        g.es["weight"] = weights
        g = g.as_undirected(combine_edges="max")

        partition = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition,
            weights="weight", resolution_parameter=resolution, seed=seed,
        )
        return np.array(partition.membership, dtype=np.int32)
    else:
        import anndata
        adata_tmp = anndata.AnnData(X=pca_emb, obsm={"X_pca": pca_emb})
        sc.pp.neighbors(adata_tmp, n_pcs=pca_emb.shape[1], n_neighbors=k,
                        random_state=seed, use_rep="X_pca")
        sc.tl.leiden(adata_tmp, resolution=resolution, random_state=seed)
        return adata_tmp.obs["leiden"].cat.codes.values.astype(np.int32)


def run_downstream_scanpy(adata, cfg: dict) -> None:
    """
    Standard scanpy downstream: PCA -> neighbors -> Leiden.
    """
    seed = cfg["seed"]
    n_comps = min(cfg["n_pcs"], adata.n_vars - 1, adata.n_obs - 1)

    if issparse(adata.X):
        if np.isnan(adata.X.data).any():
            adata.X.data = np.nan_to_num(adata.X.data, nan=0.0)
    else:
        if np.isnan(adata.X).any():
            adata.X = np.nan_to_num(adata.X, nan=0.0)

    sc.pp.pca(adata, n_comps=n_comps, random_state=seed)
    sc.pp.neighbors(adata, n_pcs=n_comps, n_neighbors=cfg["k_neighbors"],
                    random_state=seed)
    sc.tl.leiden(adata, resolution=cfg["resolution"], random_state=seed)


# ============================================================================
# PHASE 2a: CLUSTERING STABILITY
# ============================================================================

def phase2a_clustering_stability(
    adata, layer_name: str, gene_mask: np.ndarray, cfg: dict,
    *, n_bootstrap: int | None = None, seed: int | None = None,
    pca_embedding: np.ndarray | None = None,
) -> tuple[float, list[float]]:
    """
    Clustering reproducibility via bootstrap subsampling.

    PCA is computed ONCE. Each bootstrap subsamples rows from the
    precomputed embedding and runs only kNN + Leiden.

    Returns (median_ARI, list_of_pairwise_ARIs).
    """
    n_boot = n_bootstrap if n_bootstrap is not None else cfg["n_bootstrap"]
    base_seed = seed if seed is not None else cfg["seed"]

    n_cells = adata.n_obs

    if n_cells > 50_000:
        n_sub = min(int(n_cells * 0.5),50_000)
    else:
        n_sub = int(n_cells * cfg["subsample_frac"])

    if pca_embedding is None:
        pca_embedding = compute_pca_embedding(
            adata, layer_name, gene_mask, cfg["n_pcs"], seed=cfg["seed"]
        )

    all_labels = np.full((n_cells, n_boot), -1, dtype=np.int32)
    max_workers = min(cfg.get("n_jobs", 4), os.cpu_count())

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _single_bootstrap,
                b, pca_embedding, n_cells, n_sub, cfg, base_seed
            )
            for b in range(n_boot)
        ]
        for future in tqdm(as_completed(futures), total=n_boot,
                           desc="  bootstrap", leave=False):
            b, cell_idx, labels = future.result()
            if cell_idx is not None:
                all_labels[cell_idx, b] = labels

    valid = [b for b in range(n_boot) if np.any(all_labels[:, b] != -1)]
    ari_scores: list[float] = []
    for i, j in combinations(valid, 2):
        shared = (all_labels[:, i] != -1) & (all_labels[:, j] != -1)
        if shared.sum() < 50:
            continue
        ari_scores.append(
            adjusted_rand_score(all_labels[shared, i], all_labels[shared, j])
        )

    if not ari_scores:
        return np.nan, []
    return float(np.median(ari_scores)), ari_scores


# ============================================================================
# PHASE 3: HVG vs. RANDOM GENE CONTROL
# ============================================================================

def phase3_hvg_vs_random(
    adata, layer_name: str, fraction: float, hvg_mask: np.ndarray,
    cfg: dict, *, hvg_pca_embedding: np.ndarray | None = None,
    phase3_ckpt_path: Path | None = None,
) -> tuple[float, list[float]]:
    """
    Compare HVG clustering stability against random gene sets.

    Returns (hvg_stability_ARI, list_of_random_stability_ARIs).
    """
    boot_per_eval = cfg.get("n_bootstrap_phase3", 100)

    # hvg_stab, _ = phase2a_clustering_stability(
    #     adata, layer_name, hvg_mask, cfg,
    #     n_bootstrap=boot_per_eval, pca_embedding=hvg_pca_embedding,
    # )
    # Load existing Phase 3 sub-checkpoint
    p3_state = {}
    if phase3_ckpt_path is not None and phase3_ckpt_path.exists():
        try:
            with open(phase3_ckpt_path) as f:
                p3_state = json.load(f)
            log(f"      Phase 3 sub-checkpoint: {list(p3_state.keys())}")
        except Exception:
            p3_state = {}

    def _save_p3():
        if phase3_ckpt_path is not None:
            phase3_ckpt_path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp = tempfile.mkstemp(
                dir=str(phase3_ckpt_path.parent), suffix=".tmp"
            )
            try:
                with os.fdopen(fd, "w") as f:
                    json.dump(_sanitize_for_json(p3_state), f)
                os.replace(tmp, str(phase3_ckpt_path))
            except Exception:
                try:
                    os.close(fd)
                except OSError:
                    pass
                if os.path.exists(tmp):
                    os.unlink(tmp)

    # HVG stability
    if "hvg_stab" in p3_state:
        hvg_stab = p3_state["hvg_stab"]
        log(f"      HVG stability loaded from checkpoint: {hvg_stab:.4f}")
    else:
        hvg_stab, _ = phase2a_clustering_stability(
            adata, layer_name, hvg_mask, cfg,
            n_bootstrap=boot_per_eval, pca_embedding=hvg_pca_embedding,
        )
        p3_state["hvg_stab"] = hvg_stab
        _save_p3()
        log(f"      HVG stability: {hvg_stab:.4f} (checkpointed)")

    rng_masks = np.random.default_rng(cfg["seed"] + 9999)
    rand_stabs: list[float] = p3_state.get("rand_stabs", [])
    n_done = len(rand_stabs)

    if adata.n_obs > 50_000:
        cell_idx = rng_masks.choice(adata.n_obs, size=50_000, replace=False)
        adata_compute = adata[cell_idx].copy()
    else:
        adata_compute = adata.copy()

    # Advance RNG state to match where we left off
    for r in range(n_done):
        _ = random_gene_mask(adata_compute.n_vars, fraction, rng_masks)

    for r in range(n_done, cfg["n_random_repeats"]):
        rmask = random_gene_mask(adata_compute.n_vars, fraction, rng_masks)
        try:
            rand_pca = compute_pca_embedding(
                adata_compute, layer_name, rmask, cfg["n_pcs"],
                seed=cfg["seed"] + r,
            )
            rs, _ = phase2a_clustering_stability(
                adata_compute, layer_name, rmask, cfg,
                n_bootstrap=boot_per_eval, seed=cfg["seed"] + r,
                pca_embedding=rand_pca,
            )
            rand_stabs.append(rs)
            del rand_pca
        except Exception:
            rand_stabs.append(np.nan)

        p3_state["rand_stabs"] = rand_stabs
        _save_p3()
        log(f"      Random repeat {r+1}/{cfg['n_random_repeats']}: "
            f"{rand_stabs[-1]:.4f} (checkpointed)")

    return hvg_stab, rand_stabs


# ============================================================================
# PHASE 4a: SPATIAL COHERENCE
# ============================================================================

def phase4a_spatial_coherence(
    labels: np.ndarray, spatial_coords: np.ndarray, k: int = 15,
) -> float:
    """
    Fraction of each cell's k spatial nearest neighbours sharing its
    cluster label, averaged over all cells. Fully vectorized.
    """
    tree = cKDTree(spatial_coords)
    _, nn_idx = tree.query(spatial_coords, k=k + 1, workers=-1)
    nn_idx = nn_idx[:, 1:]

    nn_labels = labels[nn_idx]
    matches = nn_labels == labels[:, None]
    per_cell = matches.sum(axis=1) / k
    return float(per_cell.mean())


# ============================================================================
# PHASE 4b: MARKER GENE DISCRIMINABILITY (against external reference labels)
# ============================================================================

def load_reference_labels(reference_path: str, obs_names) -> pd.Series | None:
    """
    Load external cell-type labels from an InSituType annotation CSV.

    The reference annotation is typically the InSituType output for
    "none" normalization (raw counts) at 100% HVG, providing a
    pipeline-independent set of cell-type labels.

    Parameters
    ----------
    reference_path : str
        Path to the InSituType annotation CSV. Expected to have an
        index column matching cell IDs and a 'sup.clust' column
        with cell-type assignments.
    obs_names : pd.Index
        Cell IDs from the AnnData object, used for alignment.

    Returns
    -------
    pd.Series or None
        Cell-type labels aligned to obs_names. Cells not present in
        the reference are assigned NaN. Returns None if loading fails.
    """
    if not reference_path or not os.path.exists(reference_path):
        return None

    try:
        ref_df = pd.read_csv(reference_path, index_col=0)
    except Exception as e:
        log(f"  WARNING: failed to load reference annotations: {e}")
        return None

    # Identify the label column
    if 'sup.clust' in ref_df.columns:
        label_col = 'sup.clust'
    elif len(ref_df.columns) == 1:
        label_col = ref_df.columns[0]
    else:
        log(f"  WARNING: reference CSV has unexpected columns: {ref_df.columns.tolist()}")
        return None

    ref_labels = ref_df[label_col]
    ref_labels.index = ref_labels.index.astype(str)

    # Align to adata obs_names
    aligned = ref_labels.reindex(obs_names.astype(str))
    n_matched = aligned.notna().sum()
    n_total = len(aligned)
    log(f"  Reference labels loaded: {n_matched}/{n_total} cells matched")

    if n_matched < 50:
        log(f"  WARNING: too few matched cells ({n_matched}), skipping reference AUROC")
        return None

    return aligned


def phase4b_marker_auroc(adata_sub, layer_name: str,
                         reference_labels: pd.Series,
                         n_top_markers: int = 5):
    """
    Compute marker-gene discriminability against external reference labels.

    For each reference cell type, finds the top marker genes by Wilcoxon
    rank-sum test, then computes one-vs-rest AUROC for each marker using
    the pipeline's normalized expression.

    Parameters
    ----------
    adata_sub : anndata.AnnData
        Subsampled AnnData with the pipeline's normalized expression in .X.
        Must have the same obs index as the reference_labels.
    layer_name : str
        Name of the normalization layer (used for dense extraction).
    reference_labels : pd.Series
        External cell-type labels aligned to adata_sub.obs_names.
        Cells with NaN labels are excluded.
    n_top_markers : int
        Number of top DE genes per cell type to evaluate (default 5).

    Returns
    -------
    tuple of (float, list[float])
        (mean_auroc, list_of_individual_aurocs).
        Returns (np.nan, []) if computation fails.
    """
    adata_w = adata_sub.copy()

    # Assign reference labels and filter to cells with valid labels
    aligned_labels = reference_labels.reindex(adata_w.obs_names.astype(str))
    valid_mask = aligned_labels.notna()
    if valid_mask.sum() < 50:
        return np.nan, []

    adata_w = adata_w[valid_mask.values].copy()
    adata_w.obs["ref_celltype"] = aligned_labels[valid_mask].values.astype(str)

    if layer_name != "X":
        adata_w.X = get_dense_layer(adata_w, layer_name)

    unique_labels = np.unique(adata_w.obs["ref_celltype"])
    if len(unique_labels) < 2:
        return np.nan, []

    # Wilcoxon DE against reference cell types
    try:
        sc.tl.rank_genes_groups(
            adata_w, groupby="ref_celltype", method="wilcoxon",
            n_genes=min(20, adata_w.n_vars),
        )
    except Exception:
        return np.nan, []

    labels = adata_w.obs["ref_celltype"].values.astype(str)
    X_dense = get_dense_layer(adata_w, "X")
    gene_list = list(adata_w.var_names)
    result = adata_w.uns["rank_genes_groups"]

    aurocs: list[float] = []
    for group in unique_labels:
        try:
            top_genes = result["names"][group][:n_top_markers]
        except (KeyError, IndexError):
            continue
        binary = (labels == group).astype(np.int32)
        for gene in top_genes:
            if gene not in gene_list:
                continue
            expr = X_dense[:, gene_list.index(gene)]
            if expr.std() == 0:
                continue
            try:
                aurocs.append(roc_auc_score(binary, expr))
            except ValueError:
                continue

    if not aurocs:
        return np.nan, []
    return float(np.mean(aurocs)), aurocs


# ============================================================================
# ATOMIC CHECKPOINT I/O
# ============================================================================

def save_checkpoint(results: list[dict], path: Path) -> None:
    """Write results atomically (write-to-temp, then rename)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp",
                               prefix=path.stem)
    try:
        pd.DataFrame(results).to_csv(tmp, index=False)
        os.close(fd)
        os.replace(tmp, str(path))
    except Exception:
        os.close(fd)
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def load_checkpoint(path: Path) -> pd.DataFrame | None:
    """Load a checkpoint CSV. Returns None if absent or corrupt."""
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


# ============================================================================
# PHASE-LEVEL CHECKPOINT I/O
# Saves partial results within a fraction so that phases (2a, 3, 4a, 4b)
# survive 6-hour wall-time kills without re-running completed phases.
# ============================================================================
 
def _phase_ckpt_path(ckpt_base: Path, layer_safe: str,
                     frac: float) -> Path:
    """Path to the per-fraction phase checkpoint JSON."""
    return (ckpt_base / "phase_checkpoints"
            / f"{layer_safe}_frac{frac:.2f}.json")
 
 
def _pca_cache_path(ckpt_base: Path, layer_safe: str,
                    frac: float) -> Path:
    """Path to the cached PCA embedding (.npy)."""
    return (ckpt_base / "phase_checkpoints"
            / f"{layer_safe}_frac{frac:.2f}_pca.npy")

 
def _sanitize_for_json(obj):
    """Convert numpy types to Python natives for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj
 
 
def save_phase_ckpt(path: Path, phase_data: dict) -> None:
    """Save phase-level checkpoint as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(_sanitize_for_json(phase_data), f)
        os.replace(tmp, str(path))
    except Exception:
        try:
            os.close(fd)
        except OSError:
            pass
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise

 
def load_phase_ckpt(path: Path) -> dict:
    """Load phase-level checkpoint. Returns empty dict if absent."""
    if not path.exists():
        return {}
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return {}
 
 
def save_pca_cache(path: Path, embedding: np.ndarray) -> None:
    """Cache PCA embedding to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(str(path), embedding)
 
 
def load_pca_cache(path: Path) -> np.ndarray | None:
    """Load cached PCA embedding. Returns None if absent."""
    if not path.exists():
        return None
    try:
        return np.load(str(path))
    except Exception:
        return None


# ============================================================================
# MAIN: PROCESS ONE (DATASET, LAYER) PAIR
# ============================================================================

def run_single_layer(
    zarr_path: str, nextflow_path: str, technology: str,
    dataset_name: str, layer_name: str, cfg: dict,
) -> None:
    """
    Load one dataset, benchmark one normalization layer across all
    HVG fractions, checkpoint after each fraction.
    """
    log(f"{'='*60}")
    log(f"Dataset: {dataset_name}  |  Layer: {layer_name}")
    log(f"{'='*60}")

    # -- Load data --------------------------------------------------------
    sd_obj = load_spatialdata(zarr_path)
    adata = sd_obj.tables["filtered"].copy()

    try:
        spatial_coords = get_spatial_coordinates(adata, sd_obj)
        has_spatial = True
        log(f"  Spatial coords: {spatial_coords.shape}")
    except ValueError as e:
        log(f"  WARNING: {e} -- spatial coherence will be skipped")
        has_spatial = False
        spatial_coords = None

    del sd_obj
    gc.collect()

    norm_map = discover_norm_layers(nextflow_path)
    load_layer(adata, layer_name, norm_map)
    norm, trans = _parse_layer(layer_name)

    # -- Load external reference labels for Phase 4b ----------------------
    ref_labels = None
    if cfg.get("reference_annotation"):
        log(f"  Loading reference labels: {cfg['reference_annotation']}")
        ref_labels = load_reference_labels(
            cfg["reference_annotation"], adata.obs_names
        )
        if ref_labels is None:
            log("  WARNING: reference labels unavailable, Phase 4b will be skipped")
    else:
        log("  No --reference_annotation provided, Phase 4b will be skipped")

    # -- Checkpoint path --------------------------------------------------
    # If checkpoint_dir is provided, write checkpoints there (persistent
    # across Nextflow work directory changes). Otherwise fall back to
    # output_dir (inside work directory, lost on re-run).
    layer_safe = layer_name.replace("/", "_")
    frac_tag = "_".join(f"{f:.2f}" for f in cfg["hvg_fractions"])
    ckpt_base = cfg["checkpoint_dir"] if cfg.get("checkpoint_dir") else cfg["output_dir"]
    ckpt_path = (
        Path(ckpt_base) / "checkpoints"
        / f"{layer_safe}_{frac_tag}.csv"
    )
    log(f"  Checkpoint path: {ckpt_path}")

    # Layer-level tmp parent. Each fraction gets its OWN subdir below this
    # (frac{f}/), so concurrent fractions of the same layer never share a
    # directory and can never delete each other's checkpoints.
    tmp_base = (
        Path(ckpt_base) / "tmp" / layer_safe
    )
    log(f"  Phase tmp dir: {tmp_base}")

    existing_df = load_checkpoint(ckpt_path)
    completed_fracs: set[float] = set()
    results: list[dict] = []
    if existing_df is not None and len(existing_df) > 0:
        results = existing_df.to_dict("records")
        completed_fracs = set(existing_df["hvg_fraction"].unique())
        log(f"  Resuming: fractions already done = {completed_fracs}")

    # -- Main loop over HVG fractions -------------------------------------
    for frac in cfg["hvg_fractions"]:
        if frac in completed_fracs:
            log(f"  Skipping frac={frac:.0%} (checkpoint exists)")
            continue

        try:
            hvg_mask = load_hvg_mask(adata, layer_name, frac, nextflow_path)
            n_selected = int(hvg_mask.sum())
        except Exception as e:
            log(f"  ERROR selecting HVGs @ {frac}: {e}")
            results.append({
                "dataset": dataset_name, "technology": technology,
                "layer": layer_name, "normalization": norm,
                "transformation": trans, "hvg_fraction": frac,
                "error": str(e),
            })
            save_checkpoint(results, ckpt_path)
            continue

        log(f"  frac={frac:.0%} ({n_selected}/{adata.n_vars} genes)")

        # -- Per-fraction tmp dir + checkpoint paths ------------------------
        # Scoping tmp to the fraction is what prevents a sibling fraction's
        # cleanup from deleting this fraction's live checkpoints.
        frac_tmp = tmp_base / f"frac{frac:.2f}"
        ph_ckpt = _phase_ckpt_path(frac_tmp, layer_safe, frac)
        pca_npy = _pca_cache_path(frac_tmp, layer_safe, frac)
 
        phase_done = load_phase_ckpt(ph_ckpt)

        if phase_done:
            log(f"    Resuming: phases done = "
                f"{[k for k in phase_done if k != '_row']}")
 
        row: dict = phase_done.get("_row", {})
        if not row:
            row = {
                "dataset": dataset_name, "technology": technology,
                "layer": layer_name, "normalization": norm,
                "transformation": trans, "hvg_fraction": frac,
                "n_genes_selected": n_selected, "n_genes_total": adata.n_vars,
            }

        # row: dict = {
        #     "dataset": dataset_name, "technology": technology,
        #     "layer": layer_name, "normalization": norm,
        #     "transformation": trans, "hvg_fraction": frac,
        #     "n_genes_selected": n_selected, "n_genes_total": adata.n_vars,
        # }

        # ---- Precompute PCA -----------------------------------------------
        hvg_pca = None
        if "pca" in phase_done:
            hvg_pca = load_pca_cache(pca_npy)
            if hvg_pca is not None:
                log(f"    PCA loaded from cache: {hvg_pca.shape}")
            else:
                log("    PCA cache file missing, recomputing...")
 
        if hvg_pca is None:
            try:
                log("    Computing PCA ...")
                hvg_pca = compute_pca_embedding(
                    adata, layer_name, hvg_mask, cfg["n_pcs"], seed=cfg["seed"]
                )
                log(f"      shape: {hvg_pca.shape}")
                save_pca_cache(pca_npy, hvg_pca)
                phase_done["pca"] = True
                phase_done["_row"] = row
                save_phase_ckpt(ph_ckpt, phase_done)
                log("      PCA cached to disk")
            except Exception as e:
                log(f"      ERROR: {e}")

        # ---- Phase 2a: Clustering stability --------------------------------
        if "phase2a" not in phase_done:
            try:
                log("    Phase 2a: Clustering stability ...")
                stab_med, stab_all = phase2a_clustering_stability(
                    adata, layer_name, hvg_mask, cfg, pca_embedding=hvg_pca,
                )
                row["clustering_stability_median_ari"] = stab_med
                row["clustering_stability_scores"] = json.dumps(
                    [round(s, 4) for s in stab_all[:50]]
                )
                if stab_all:
                    mc_se = np.std(stab_all) / np.sqrt(len(stab_all))
                    row["clustering_stability_mc_se"] = round(mc_se, 6)
                log(f"      median ARI = {stab_med:.4f}")
 
                phase_done["phase2a"] = True
                phase_done["_row"] = row
                save_phase_ckpt(ph_ckpt, phase_done)
            except Exception as e:
                log(f"      ERROR 2a: {e}")
                row["clustering_stability_median_ari"] = np.nan
        else:
            log("    Phase 2a: skipped (checkpoint)")

        # ---- Phase 3: HVG vs random (skip for frac=1.0) -------------------
        if frac < 1.0 and "phase3" not in phase_done:
            try:
                log("    Phase 3: HVG vs. random control ...")
                p3_ckpt = frac_tmp / "phase_checkpoints" / f"{layer_safe}_frac{frac:.2f}_phase3.json"
                hvg_stab, rand_stabs = phase3_hvg_vs_random(
                    adata, layer_name, frac, hvg_mask, cfg,
                    hvg_pca_embedding=hvg_pca,
                    phase3_ckpt_path=p3_ckpt,
                )
                row["hvg_stability"] = hvg_stab
                row["random_stability_mean"] = float(np.nanmean(rand_stabs))
                row["random_stability_std"] = float(np.nanstd(rand_stabs))
                row["random_stability_values"] = json.dumps(
                    [round(s, 4) for s in rand_stabs if not np.isnan(s)]
                )
                valid_rs = [s for s in rand_stabs if not np.isnan(s)]
                if valid_rs:
                    row["random_stability_mc_se"] = float(
                        np.std(valid_rs) / np.sqrt(len(valid_rs))
                    )
                log(f"      HVG={hvg_stab:.4f}, "
                    f"Random={np.nanmean(rand_stabs):.4f} "
                    f"+/- {np.nanstd(rand_stabs):.4f}")
 
                phase_done["phase3"] = True
                phase_done["_row"] = row
                save_phase_ckpt(ph_ckpt, phase_done)

                # Clean Phase 3 sub-checkpoint
                if p3_ckpt.exists():
                    p3_ckpt.unlink()
            except Exception as e:
                log(f"      ERROR 3: {e}")
        elif frac < 1.0:
            log("    Phase 3: skipped (checkpoint)")
        
        # ---- Phase 4: Biological coherence ---------------------------------
        if "phase4" not in phase_done:
            adata_clust = None
            try:
                rng = np.random.default_rng(cfg["seed"])
                n_sub = int(adata.n_obs * cfg["subsample_frac"])
                cell_idx = np.sort(rng.choice(adata.n_obs, size=n_sub,
                                              replace=False))
 
                adata_sub = adata[cell_idx, :][:, hvg_mask].copy()
                if layer_name != "X":
                    adata_sub.X = get_dense_layer(adata_sub, layer_name)
                run_downstream_scanpy(adata_sub, cfg)
                adata_clust = adata_sub
                row["n_clusters"] = int(adata_clust.obs["leiden"].nunique())
            except Exception as e:
                log(f"      ERROR building clusters for Phase 4: {e}")
 
            # Phase 4a: Spatial coherence
            if adata_clust is not None and has_spatial:
                try:
                    log("    Phase 4a: Spatial coherence ...")
                    labels_int = adata_clust.obs["leiden"].cat.codes.values
                    coords_sub = spatial_coords[cell_idx]
                    sp_coh = phase4a_spatial_coherence(
                        labels_int, coords_sub, k=cfg["k_neighbors"],
                    )
                    row["spatial_coherence"] = sp_coh
                    log(f"      coherence = {sp_coh:.4f}")
                except Exception as e:
                    log(f"      ERROR 4a: {e}")
                    row["spatial_coherence"] = np.nan
 
            # Phase 4b: Marker gene discriminability (against reference labels)
            if adata_clust is not None and ref_labels is not None:
                try:
                    log("    Phase 4b: Marker gene discriminability (ref labels) ...")
                    auroc_mean, _ = phase4b_marker_auroc(
                        adata_clust, layer_name, ref_labels
                    )
                    row["marker_auroc_mean"] = auroc_mean
                    log(f"      mean AUROC = {auroc_mean:.4f}")
                except Exception as e:
                    log(f"      ERROR 4b: {e}")
                    row["marker_auroc_mean"] = np.nan
 
            del adata_clust
            phase_done["phase4"] = True
            phase_done["_row"] = row
            save_phase_ckpt(ph_ckpt, phase_done)
        else:
            log("    Phase 4: skipped (checkpoint)")
 
        del hvg_pca
        gc.collect()

        # ---- Checkpoint this fraction --------------------------------------
        results.append(row)
        save_checkpoint(results, ckpt_path)

        # Clean up this fraction's own tmp tree. Because frac_tmp is unique
        # to this fraction, this can never touch another fraction's files.
        if frac_tmp.exists():
            try:
                shutil.rmtree(frac_tmp)
            except OSError:
                pass
        log(f"    Fraction {frac:.0%} complete, phase tmp cleaned up")
 
    log(f"Finished: {dataset_name} / {layer_name} ({len(results)} rows)")
    log(f"Checkpoint: {ckpt_path}")

    # Also write final results to the work directory output_dir so Nextflow
    # can capture it as a process output (checkpoint_dir may be external)
    if cfg.get("checkpoint_dir") and str(ckpt_path.parent) != str(Path(cfg["output_dir"]) / "checkpoints"):
        final_path = (
            Path(cfg["output_dir"]) / "checkpoints"
            / f"{layer_safe}_{frac_tag}.csv"
        )
        save_checkpoint(results, final_path)
        log(f"Final output: {final_path}")

    # Remove the layer tmp parent ONLY if empty. Never rmtree it: a sibling
    # fraction of the same layer may be running concurrently and still owns
    # its frac{f}/ subdir below this path.
    #if tmp_base.exists():
    #    try:
    #        tmp_base.rmdir()
    #        log(f"  Cleaned up empty tmp dir: {tmp_base}")
    #    except OSError:
    #        pass  # other fractions still have subdirs here


# ============================================================================
# CLI
# ============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="BISTRO: HVG Selection Benchmarking (single dataset+layer)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--zarr", required=True,
                   help="Path to the QC'd .zarr file")
    p.add_argument("--nextflow_output", required=True,
                   help="Path to the Nextflow output directory")
    p.add_argument("--technology", required=True,
                   choices=["CosMx", "Xenium", "MERFISH"],
                   help="iST platform")
    p.add_argument("--dataset_name", default=None,
                   help="Dataset name (default: inferred from zarr path)")
    p.add_argument("--layer", type=str, default=None,
                   help="Normalization layer name. "
                        "Required unless --list_layers is used.")
    p.add_argument("--list_layers", action="store_true",
                   help="Print available layers and exit.")
    p.add_argument("--output_dir", required=True,
                   help="Directory to write output files")
    p.add_argument("--n_bootstrap", type=int, default=15)
    p.add_argument("--n_bootstrap_phase3", type=int, default=100,
                   help="Bootstrap iterations per evaluation in Phase 3 "
                        "HVG vs random comparison (default 100)")
    p.add_argument("--subsample_frac", type=float, default=0.80)
    p.add_argument("--n_random_repeats", type=int, default=5)
    p.add_argument("--k_neighbors", type=int, default=15)
    p.add_argument("--resolution", type=float, default=0.5)
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--hvg_fractions", type=float, nargs="+",
                   default=[0.10, 0.25, 0.50, 0.75, 1.00])
    p.add_argument("--n_jobs", type=int, default=4)
    p.add_argument("--checkpoint_dir", default=None,
                   help="Persistent directory for checkpoints that survives "
                        "Nextflow work directory changes. If not provided, "
                        "checkpoints are written inside --output_dir (only "
                        "useful if that path is persistent).")
    p.add_argument("--reference_annotation", default=None,
                   help="Path to InSituType annotation CSV for the reference "
                        "pipeline (e.g. 'none' normalization, 100%% HVG). "
                        "Used for marker-gene discriminability (Phase 4b). "
                        "If not provided, Phase 4b is skipped.")
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    dataset_name = args.dataset_name
    if dataset_name is None:
        dataset_name = os.path.basename(args.zarr).replace('.zarr', '')

    norm_map = discover_norm_layers(args.nextflow_output)

    if args.list_layers:
        layers = sorted(norm_map.keys())
        print(f"Dataset: {dataset_name}")
        print(f"Layers ({len(layers)}):")
        for ln in layers:
            print(f"  {ln}")
        sys.exit(0)

    if args.layer is None:
        print("ERROR: --layer is required (use --list_layers to see options).")
        sys.exit(1)

    if args.layer not in norm_map:
        print(f"ERROR: Layer '{args.layer}' not found.")
        print(f"Available: {sorted(norm_map.keys())}")
        sys.exit(1)

    cfg = {
        "output_dir":           args.output_dir,
        "checkpoint_dir":       args.checkpoint_dir,
        "reference_annotation": args.reference_annotation,
        "n_bootstrap":          args.n_bootstrap,
        "subsample_frac":       args.subsample_frac,
        "n_random_repeats":     args.n_random_repeats,
        "k_neighbors":          args.k_neighbors,
        "resolution":           args.resolution,
        "n_pcs":                args.n_pcs,
        "seed":                 args.seed,
        "hvg_fractions":        args.hvg_fractions,
        "n_jobs":               args.n_jobs,
        "n_bootstrap_phase3":   args.n_bootstrap_phase3,
    }

    log(f"Config:\n{json.dumps(cfg, indent=2)}")
    run_single_layer(
        zarr_path=args.zarr,
        nextflow_path=args.nextflow_output,
        technology=args.technology,
        dataset_name=dataset_name,
        layer_name=args.layer,
        cfg=cfg,
    )
    log("Done.")


if __name__ == "__main__":
    main()