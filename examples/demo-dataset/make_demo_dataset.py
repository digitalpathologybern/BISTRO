#!/usr/bin/env python
"""
BISTRO demo dataset generator
=============================

Builds a small, fully simulated iST dataset that exercises every stage of the
BISTRO pipeline on a laptop in minutes, and doubles as the *executable
specification* of the ``filtered.zarr`` input schema documented in the README.

The simulation plants a known ground truth so the demo verifies correctness,
not merely that the code runs:

  * a per-FOV library-size offset with known variance  -> the FOV batch effect
    that ``bin/evaluation/batch_effect.py`` estimates as tau^2 (Equation 2)
  * a systematic decline in that offset across the FOV acquisition order
    -> the acquisition drift reported as ``drift_slope`` / ``drift_pearson_r``
  * a designated set of cell-type marker genes -> the genes HVG selection
    should preferentially recover
  * five cell types with distinct expression profiles, also written out as the
    InSituType reference -> annotation should recover the true labels

The realised values of all of these are written to ``ground_truth.json`` so a
demo run can be checked against them.

Counts follow a negative-binomial model, matching the mean-variance
relationship (Var = mu + mu^2 / theta) that the transformation analysis in
Section 2.3 is designed to characterise.

Outputs (into --output_dir)
---------------------------
  demo_filtered.zarr              SpatialData archive, table 'filtered'
  demo_tissue_annotations.csv     tissue regions, one row per cell
  demo_reference_profiles.csv     genes x cell types, input to the R step
  demo_color_list.txt             colour list for create_colormap.R
  ground_truth.json               simulation parameters and realised values

``demo_reference_profiles.csv`` is converted to the ``.RData`` that InSituType
expects by ``make_demo_reference.R``; ``run_demo.sh`` runs both in order.

Usage
-----
    python make_demo_dataset.py --output_dir .
"""

import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from spatialdata import SpatialData
from spatialdata.models import TableModel


# ============================================================================
# Simulation constants
# ============================================================================

N_CELLS = 2000
N_GENES = 200
# 4 x 4 = 16 FOVs. The number of FOVs is the number of groups the mixed-effects
# model has to estimate tau^2 from, so it drives how precisely the injected
# batch effect can be recovered; 16 keeps REML shrinkage modest while leaving
# ~125 cells per FOV, ample for estimating a per-FOV mean library size.
FOV_GRID = (4, 4)
FOV_SIZE_UM = 510.72           # CosMx FOV edge length, matches get_fov_size()
N_TISSUE_REGIONS = 3
SEED = 42

# Negative-binomial dispersion. Var = mu + mu^2 / THETA.
THETA = 10.0

# Per-FOV log library-size offset: u_f = DRIFT_SLOPE * centred_index + noise.
# DRIFT_SLOPE produces a systematic signal decline across the acquisition
# order; FOV_NOISE_SD is the unstructured FOV-to-FOV component. Together they
# define the tau^2 that the batch-effect model should recover.
#
# DRIFT_SLOPE is set so the systematic component dominates the unstructured
# one: its spread across 16 FOVs is 0.035 * sd(0..15) = 0.161, against a noise
# SD of 0.150, giving an expected correlation of about -0.73 between the FOV
# offset and the acquisition index. A weaker drift makes the demonstration
# hostage to the draw -- at -0.018 the seed-42 draw came out at r = -0.16,
# in the bottom 4% of draws, and the drift analysis had almost nothing to
# find. Over 16 FOVs this is a ~41% decline in signal from first tile to last,
# which is strong but within what a poor acquisition run really shows.
DRIFT_SLOPE = -0.035
FOV_NOISE_SD = 0.150

# Residual per-cell variation in library size, on the log scale.
CELL_NOISE_SD = 0.150

CELL_TYPES = ["Epithelial", "Tcell", "Bcell", "Myeloid", "Stromal"]
N_MARKERS_PER_TYPE = 8         # 5 x 8 = 40 marker genes = the HVG ground truth
MARKER_LOG_FC = 1.6            # log-scale enrichment of a marker in its type

TISSUE_REGIONS = ["Tumor", "Stroma", "Immune"]
N_TISSUE_GENES = 20            # genes carrying a tissue-region effect
TISSUE_LOG_FC = 0.7

# Real human gene symbols, deliberately free of '-', '_', ':', '/' and spaces.
# The R stages normalise those characters to '.' (see clean_gene_names() in
# bin/annotation/run_insitutype.R), so avoiding them keeps the demo's gene
# names identical across every stage.
GENE_PANEL = [
    # Epithelial / tumour
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "CDH1", "CDX2",
    "VIL1", "MUC2", "LGR5", "OLFM4", "TFF3", "REG4", "SPINK4", "AGR2",
    "CEACAM5", "CEACAM6", "MKI67", "PCNA", "TOP2A", "CCND1", "MYC", "TP53",
    "CTNNB1", "AXIN2", "ASCL2", "SOX9", "HNF4A", "KLF5",
    # T / NK
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B",
    "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "FOXP3", "IL2RA", "CTLA4",
    "PDCD1", "LAG3", "TIGIT", "HAVCR2", "GZMA", "GZMB", "GZMK", "PRF1",
    "NKG7", "KLRD1", "KLRB1", "NCAM1", "FCGR3A", "TRDC",
    # B / plasma
    "CD19", "MS4A1", "CD79A", "CD79B", "TCL1A", "IGHM", "IGHD", "IGHG1",
    "IGHA1", "JCHAIN", "MZB1", "XBP1", "PRDM1", "IRF4", "TNFRSF17",
    # Myeloid
    "LYZ", "CD14", "CD68", "CD163", "MRC1", "MSR1", "CSF1R", "ITGAM",
    "ITGAX", "FCGR1A", "S100A8", "S100A9", "S100A12", "FCN1", "VCAN",
    "CLEC9A", "CLEC10A", "CD1C", "LAMP3", "IRF8", "BATF3", "FLT3", "SIRPA",
    # Stromal / fibroblast / muscle
    "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL4A1", "FN1",
    "LUM", "DCN", "POSTN", "THY1", "PDGFRA", "PDGFRB", "ACTA2", "TAGLN",
    "MYH11", "DES", "RGS5", "NOTCH3", "CSPG4",
    # Endothelial / lymphatic
    "PECAM1", "VWF", "CDH5", "CLDN5", "PLVAP", "ESM1", "ANGPT2", "KDR",
    "FLT1", "TEK", "LYVE1", "PROX1", "PDPN", "CCL21",
    # Cytokines / chemokines
    "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL8", "CCL2",
    "CCL5", "CCL19", "CCL20", "IL1B", "IL6", "IL10", "IL17A", "IFNG",
    "TNF", "TGFB1", "TGFB2", "VEGFA", "HIF1A",
    # Antigen presentation
    "CD74", "CIITA", "TAP1", "TAP2", "B2M", "PSMB9",
    # Signalling / oncogenes
    "STAT1", "STAT3", "JAK1", "JAK2", "NFKB1", "RELA", "MAPK1", "AKT1",
    "PIK3CA", "PTEN", "MTOR", "EGFR", "ERBB2", "ERBB3", "MET", "KRAS",
    "BRAF", "SMAD4", "TGFBR2", "APC",
    # Apoptosis / cell cycle
    "BCL2", "BAX", "CASP3", "CASP8", "FAS", "FASLG", "TNFSF10", "BIRC5",
    "MCL1", "CDKN1A", "CDKN2A", "RB1", "E2F1", "CCNB1", "CDK1", "AURKA",
    "PLK1", "BUB1",
    # Matrix remodelling
    "SPP1", "MMP1", "MMP2", "MMP7", "MMP9", "TIMP1", "TIMP2", "SERPINE1",
    "PLAU", "PLAUR",
]


# ============================================================================
# Simulation
# ============================================================================

def build_gene_panel(n_genes):
    """Return the first ``n_genes`` symbols of the curated panel."""
    if len(GENE_PANEL) < n_genes:
        raise ValueError(
            f"Gene panel holds {len(GENE_PANEL)} symbols, {n_genes} requested."
        )
    genes = GENE_PANEL[:n_genes]
    if len(set(genes)) != len(genes):
        dupes = sorted({g for g in genes if genes.count(g) > 1})
        raise ValueError(f"Duplicate gene symbols in panel: {dupes}")
    return genes


def assign_cells(rng, n_cells, grid, fov_size_um):
    """
    Place cells on a regular FOV grid and assign cell types and tissue regions.

    FOVs are numbered in row-major (raster) order, matching the acquisition
    order that the drift analysis regresses against. Tissue regions are laid
    out as horizontal bands so that they are spatially contiguous, which is
    what the spatial-coherence phase of the HVG benchmark measures.

    Returns
    -------
    dict of np.ndarray
        fov, x_local_um, y_local_um, x_global_um, y_global_um,
        cell_type_idx, tissue_idx
    """
    n_rows, n_cols = grid
    n_fovs = n_rows * n_cols

    # Even split of cells across FOVs, remainder spread over the first FOVs.
    per_fov = np.full(n_fovs, n_cells // n_fovs, dtype=int)
    per_fov[: n_cells % n_fovs] += 1
    fov = np.repeat(np.arange(n_fovs), per_fov)

    # Local coordinates within each FOV, then global by grid offset.
    x_local = rng.uniform(0.0, fov_size_um, size=n_cells)
    y_local = rng.uniform(0.0, fov_size_um, size=n_cells)
    fov_row, fov_col = fov // n_cols, fov % n_cols
    x_global = fov_col * fov_size_um + x_local
    y_global = fov_row * fov_size_um + y_local

    # Tissue regions as concentric rings about the slide centre: a tumour core,
    # a stromal rim, and an immune periphery.
    #
    # Rings rather than bands, deliberately. The batch-effect model is
    # `log_LS ~ C(tissue_annotations) + (1|FOV)`, and the drift analysis
    # regresses the per-FOV random intercepts on the FOV index. Horizontal
    # bands would be almost collinear with a row-major FOV index (measured at
    # r = -0.84), so the tissue fixed effect would absorb the injected
    # acquisition drift and leave nothing for the drift analysis to find.
    # Radial distance is symmetric about the middle rows, so it is nearly
    # uncorrelated with acquisition order while staying spatially contiguous,
    # which the HVG spatial-coherence phase needs.
    centre_x = n_cols * fov_size_um / 2.0
    centre_y = n_rows * fov_size_um / 2.0
    radius = np.hypot(x_global - centre_x, y_global - centre_y)

    # Equal-count rings, so no region is too small to model.
    edges = np.quantile(radius, np.linspace(0, 1, N_TISSUE_REGIONS + 1)[1:-1])
    band = np.digitize(radius, edges)

    # Blur the boundaries slightly so regions are contiguous but not perfectly
    # separable, as in real tissue.
    flip = rng.random(n_cells) < 0.05
    tissue_idx = np.where(flip, rng.integers(0, N_TISSUE_REGIONS, n_cells), band)

    # Cell types correlate with tissue region but are not determined by it.
    # Rows are tissue regions (Tumor, Stroma, Immune), columns are cell types.
    type_probs = np.array([
        [0.60, 0.10, 0.05, 0.10, 0.15],   # Tumor
        [0.10, 0.10, 0.05, 0.15, 0.60],   # Stroma
        [0.05, 0.40, 0.20, 0.30, 0.05],   # Immune
    ])
    cell_type_idx = np.array([
        rng.choice(len(CELL_TYPES), p=type_probs[t]) for t in tissue_idx
    ])

    return {
        "fov": fov,
        "x_local_um": x_local,
        "y_local_um": y_local,
        "x_global_um": x_global,
        "y_global_um": y_global,
        "cell_type_idx": cell_type_idx,
        "tissue_idx": tissue_idx,
    }


def build_expression_profiles(rng, genes):
    """
    Build the per-cell-type and per-tissue log-expression profiles.

    A contiguous, disjoint block of ``N_MARKERS_PER_TYPE`` genes is assigned to
    each cell type and up-weighted by ``MARKER_LOG_FC`` in that type only.
    Those marker genes carry essentially all of the between-cell-type variance
    and constitute the HVG ground truth.

    Returns
    -------
    tuple
        (beta, gamma, markers) where beta is (n_types, n_genes) log-expression
        per cell type, gamma is (n_tissues, n_genes) additive tissue effects,
        and markers maps each cell type to its marker gene symbols.
    """
    n_genes = len(genes)
    n_types = len(CELL_TYPES)

    # Baseline log-expression. Centred near 1 count with a wide spread, which
    # reproduces the ~40-50% zero fraction and the long-tailed abundance
    # distribution characteristic of targeted iST panels.
    baseline = rng.normal(loc=np.log(1.0), scale=1.1, size=n_genes)
    beta = np.tile(baseline, (n_types, 1))

    # Disjoint marker blocks, drawn from a shuffled gene order so that marker
    # identity is not confounded with position in the panel.
    order = rng.permutation(n_genes)
    markers = {}
    for k in range(n_types):
        idx = order[k * N_MARKERS_PER_TYPE:(k + 1) * N_MARKERS_PER_TYPE]
        beta[k, idx] += MARKER_LOG_FC
        markers[CELL_TYPES[k]] = [genes[i] for i in sorted(idx)]

    # Tissue-region effects on a separate, non-overlapping set of genes.
    n_marker_total = n_types * N_MARKERS_PER_TYPE
    tissue_idx_genes = order[n_marker_total:n_marker_total + N_TISSUE_GENES]
    gamma = np.zeros((N_TISSUE_REGIONS, n_genes))
    for t in range(N_TISSUE_REGIONS):
        gamma[t, tissue_idx_genes] = rng.normal(0.0, TISSUE_LOG_FC,
                                                size=N_TISSUE_GENES)

    tissue_genes = [genes[i] for i in sorted(tissue_idx_genes)]
    return beta, gamma, markers, tissue_genes


def simulate_counts(rng, cells, beta, gamma, n_fovs):
    """
    Draw negative-binomial counts with an injected per-FOV library-size offset.

        log s_i = u_{f(i)} + e_i
        u_f     = DRIFT_SLOPE * centred_fov_index + N(0, FOV_NOISE_SD^2)
        mu_ig   = s_i * exp(beta[k(i), g] + gamma[t(i), g])
        y_ig    ~ NB(mean = mu_ig, dispersion = THETA)

    ``u_f`` is the quantity the mixed-effects model recovers as its per-FOV
    random intercept; its variance is the tau^2 of Equation 2.

    Returns
    -------
    tuple
        (counts, size_factors, u_fov)
    """
    n_cells = len(cells["fov"])

    centred = np.arange(n_fovs) - (n_fovs - 1) / 2.0
    u_fov = DRIFT_SLOPE * centred + rng.normal(0.0, FOV_NOISE_SD, size=n_fovs)

    log_s = u_fov[cells["fov"]] + rng.normal(0.0, CELL_NOISE_SD, size=n_cells)
    size_factors = np.exp(log_s)

    log_mu = (beta[cells["cell_type_idx"], :]
              + gamma[cells["tissue_idx"], :]
              + log_s[:, None])
    mu = np.exp(log_mu)

    # numpy parameterises NB by (n successes, p); mean = n(1-p)/p gives
    # p = theta / (theta + mu) for dispersion theta.
    p = THETA / (THETA + mu)
    counts = rng.negative_binomial(THETA, p).astype(np.float32)

    return counts, size_factors, u_fov


def build_obs(cells, counts, rng, genes):
    """
    Assemble the cell metadata table, i.e. ``adata.obs``.

    Every column required or consumed by the pipeline is set here; this
    function is the practical reference for the obs schema documented in the
    README.
    """
    n_cells = counts.shape[0]
    cell_ids = [f"cell_{i:05d}" for i in range(n_cells)]

    obs = pd.DataFrame(index=pd.Index(cell_ids, name=None))

    # --- FOV identifier (native, as on CosMx) --------------------------------
    obs["fov"] = cells["fov"].astype(int)

    # --- Spatial coordinates, micrometres ------------------------------------
    obs["x_local_um"] = cells["x_local_um"]
    obs["y_local_um"] = cells["y_local_um"]
    obs["x_global_um"] = cells["x_global_um"]
    obs["y_global_um"] = cells["y_global_um"]

    # --- Cell area, square micrometres (drives areaNorm) ---------------------
    obs["area_um2"] = rng.lognormal(mean=np.log(50.0), sigma=0.30, size=n_cells)

    # --- Negative-probe totals (InSituType background estimate) --------------
    # run_insitutype.R uses total_counts_Negative / 20 as the per-cell mean
    # negative-probe count, so this must be the SUM over the negative probes.
    obs["total_counts_Negative"] = rng.poisson(5.0, size=n_cells).astype(float)

    # --- Convenience QC columns, carried through but not required ------------
    obs["total_counts"] = counts.sum(axis=1)
    obs["n_genes_by_counts"] = (counts > 0).sum(axis=1)

    # --- Ground-truth labels, for demo verification only ---------------------
    # The pipeline never reads these; they let the demo check InSituType's
    # calls against the simulated truth.
    obs["true_cell_type"] = pd.Categorical(
        [CELL_TYPES[i] for i in cells["cell_type_idx"]]
    )
    obs["true_tissue_region"] = pd.Categorical(
        [TISSUE_REGIONS[i] for i in cells["tissue_idx"]]
    )

    return obs, cell_ids


def build_reference_profiles(beta, genes):
    """
    Build the InSituType reference: mean expression per gene per cell type.

    Returned on the linear scale, genes x cell types, which is the orientation
    ``run_insitutype.R`` expects from ``profile_matrix``.
    """
    return pd.DataFrame(
        np.exp(beta).T, index=genes, columns=CELL_TYPES,
    )


# ============================================================================
# Writing
# ============================================================================

def write_outputs(output_dir, adata, obs, tissue_labels, profiles,
                  ground_truth):
    """Write the zarr archive and every companion file the pipeline needs."""
    os.makedirs(output_dir, exist_ok=True)

    # --- SpatialData zarr ----------------------------------------------------
    # The table carries no region annotation: BISTRO reads expression and obs
    # only, and never dereferences shapes or images.
    table = TableModel.parse(adata)
    sdata = SpatialData(tables={"filtered": table})

    zarr_path = os.path.join(output_dir, "demo_filtered.zarr")
    if os.path.exists(zarr_path):
        import shutil
        shutil.rmtree(zarr_path)
    sdata.write(zarr_path)
    print(f"  wrote {zarr_path}")

    # --- Tissue annotations --------------------------------------------------
    # Written with the cell ID as an unnamed index column, which
    # load_tissue_annotations() picks up as 'Unnamed: 0' and uses as the merge
    # key against adata.obs.index.
    tissue_path = os.path.join(output_dir, "demo_tissue_annotations.csv")
    pd.DataFrame({"tissue_annotations": tissue_labels},
                 index=obs.index).to_csv(tissue_path)
    print(f"  wrote {tissue_path}")

    # --- InSituType reference profiles (converted to .RData by the R step) ---
    profiles_path = os.path.join(output_dir, "demo_reference_profiles.csv")
    profiles.to_csv(profiles_path)
    print(f"  wrote {profiles_path}")

    # --- Colour list for create_colormap.R -----------------------------------
    # Read with sep='\t' and header=FALSE; only the first column is used.
    colors_path = os.path.join(output_dir, "demo_color_list.txt")
    palette = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
               "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"]
    with open(colors_path, "w") as fh:
        fh.write("\n".join(palette) + "\n")
    print(f"  wrote {colors_path}")

    # --- Ground truth --------------------------------------------------------
    gt_path = os.path.join(output_dir, "ground_truth.json")
    with open(gt_path, "w") as fh:
        json.dump(ground_truth, fh, indent=2)
    print(f"  wrote {gt_path}")


# ============================================================================
# Entry point
# ============================================================================

def make_demo_dataset(output_dir, seed=SEED):
    """Generate the full demo input bundle. Deterministic given ``seed``."""
    rng = np.random.default_rng(seed)
    n_fovs = FOV_GRID[0] * FOV_GRID[1]

    print(f"BISTRO demo dataset  (seed={seed})")
    print(f"  {N_CELLS} cells x {N_GENES} genes, {n_fovs} FOVs, "
          f"{len(CELL_TYPES)} cell types, {N_TISSUE_REGIONS} tissue regions")

    genes = build_gene_panel(N_GENES)
    cells = assign_cells(rng, N_CELLS, FOV_GRID, FOV_SIZE_UM)
    beta, gamma, markers, tissue_genes = build_expression_profiles(rng, genes)
    counts, size_factors, u_fov = simulate_counts(rng, cells, beta, gamma,
                                                  n_fovs)
    obs, cell_ids = build_obs(cells, counts, rng, genes)

    var = pd.DataFrame(index=pd.Index(genes, name=None))
    adata = ad.AnnData(X=csr_matrix(counts), obs=obs, var=var)

    profiles = build_reference_profiles(beta, genes)
    tissue_labels = [TISSUE_REGIONS[i] for i in cells["tissue_idx"]]

    # ---- Realised ground truth ---------------------------------------------
    # The nominal parameters describe the generative model; the realised values
    # describe this particular draw, and are what a demo run should recover.
    lib_size = counts.sum(axis=1)
    log_lib = np.log1p(lib_size)
    fov_mean_log_lib = np.array(
        [log_lib[cells["fov"] == f].mean() for f in range(n_fovs)]
    )

    # Drift actually present in this draw. This, not the nominal DRIFT_SLOPE,
    # is what the pipeline's drift regression can recover: with 16 FOVs the
    # realised slope scatters appreciably around the nominal one.
    fov_index = np.arange(n_fovs)
    realised_drift_slope = float(np.polyfit(fov_index, u_fov, 1)[0])
    realised_drift_r = float(np.corrcoef(fov_index, u_fov)[0, 1])

    ground_truth = {
        "seed": seed,
        "n_cells": int(N_CELLS),
        "n_genes": int(N_GENES),
        "n_fovs": int(n_fovs),
        "technology_profile": "CosMx (native FOVs, 510.72 um tiles)",
        "nb_dispersion_theta": THETA,
        "batch_effect": {
            "description": (
                "Per-FOV offset u_f applied to log library size. The MELM "
                "random-intercept variance reported in "
                "*_batch_effect_summary.csv (var_random_intercept_reml) "
                "estimates var(u_f) for the 'none' normalization."
            ),
            "drift_slope_per_fov": DRIFT_SLOPE,
            "fov_noise_sd": FOV_NOISE_SD,
            "cell_noise_sd": CELL_NOISE_SD,
            "u_fov": [round(float(v), 6) for v in u_fov],
            "realised_var_u_fov": round(float(np.var(u_fov, ddof=1)), 6),
            "realised_drift_slope": round(realised_drift_slope, 6),
            "realised_drift_pearson_r": round(realised_drift_r, 4),
            "realised_fov_mean_log1p_library_size": [
                round(float(v), 6) for v in fov_mean_log_lib
            ],
        },
        "hvg": {
            "description": (
                "Marker genes carry the between-cell-type variance and are "
                "the genes HVG selection should preferentially recover."
            ),
            "n_marker_genes": int(len(CELL_TYPES) * N_MARKERS_PER_TYPE),
            "marker_log_fold_change": MARKER_LOG_FC,
            "markers_by_cell_type": markers,
            "all_marker_genes": sorted(
                {g for gs in markers.values() for g in gs}
            ),
            "tissue_effect_genes": tissue_genes,
        },
        "cell_types": {
            "labels": CELL_TYPES,
            "counts": {
                ct: int((obs["true_cell_type"] == ct).sum())
                for ct in CELL_TYPES
            },
            "note": (
                "obs['true_cell_type'] holds the simulated label. Compare "
                "against annotation/1.0_HVG_none_annotation.csv; labels are "
                "recovered up to a permutation, so score with ARI."
            ),
        },
        "tissue_regions": {
            "labels": TISSUE_REGIONS,
            "counts": {
                tr: int((obs["true_tissue_region"] == tr).sum())
                for tr in TISSUE_REGIONS
            },
        },
        "library_size": {
            "mean": round(float(lib_size.mean()), 2),
            "min": int(lib_size.min()),
            "max": int(lib_size.max()),
        },
    }

    print(f"  library size: mean {lib_size.mean():.0f}, "
          f"range [{lib_size.min()}, {lib_size.max()}]")
    print(f"  injected var(u_fov) = "
          f"{ground_truth['batch_effect']['realised_var_u_fov']:.5f}")
    print(f"  sparsity: {100 * (counts == 0).mean():.1f}% zeros")

    write_outputs(output_dir, adata, obs, tissue_labels, profiles,
                  ground_truth)
    print("Done.")
    return ground_truth


def build_parser():
    p = argparse.ArgumentParser(
        description="Generate the BISTRO demo dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--output_dir", default=".",
                   help="Directory to write the demo input bundle into")
    p.add_argument("--seed", type=int, default=SEED,
                   help="RNG seed; the output is deterministic given this")
    return p


def main():
    args = build_parser().parse_args()
    make_demo_dataset(args.output_dir, seed=args.seed)


if __name__ == "__main__":
    main()
