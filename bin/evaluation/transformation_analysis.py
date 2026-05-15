#!/usr/bin/env python
"""
BISTRO Transformation Analysis
===============================

Evaluates variance stabilization across normalization-transformation
combinations by computing per-gene mean/variance/PCA statistics and
aggregating them into metrics described in Section 2.3 of the manuscript.

Two phases:
  Phase 1 (per-layer):
      Load each normalized expression matrix, apply transformations
      (sqrt, log(x+c) with multiple pseudocounts, identity), compute
      per-gene mean, variance, and PCA loadings (5 components) for
      every derived normalization-transformation combination.

  Phase 2 (aggregation):
      Read the per-gene statistics and compute four variance stabilization
      metrics per dataset x normalization x transformation:
        1. Log-log slope: log10(var) ~ log10(mean)
        2. Spearman rho: rank correlation between mean and variance
        3. LOWESS R-squared: smooth mean-variance dependence
        4. Binned CV: coefficient of variation of median variance across
           equal-count mean-expression bins (Equation 7)

      Plus the PC1-mean correlation (Equation 6).

Outputs
-------
  {dataset}_mean_variance_stats.csv
      Per-gene statistics for every normalization-transformation combination.
      Columns: geneID, mean, variance, PC1-PC5, layer, dataset, technology.

  {dataset}_transformation_summary.csv
      Aggregated metrics per normalization-transformation combination.
      Columns: normalization, transformation, dataset, slope, spearman,
      lowess_r2, binned_cv.

  {dataset}_pc1_correlation.csv
      PC1-mean correlation per normalization-transformation combination.
      Columns: normalization, transformation, dataset, correlation.

Usage
-----
    python transformation_analysis.py \\
        --zarr /path/to/dataset.zarr \\
        --nextflow_output /path/to/nextflow_output/ \\
        --output_dir /path/to/output/ \\
        --pseudocounts 0.01 0.1 0.5 1 10
"""

import argparse
import gc
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from scipy.stats import linregress, spearmanr

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'utils'))
from helpers import parse_layer, RAW_ONLY_NORMS, COMPOUND_NORMS
from data_loader import (
    load_spatialdata, discover_norm_layers, load_layer_into_spatialdata,
    unload_layer,
)


# ============================================================================
# Transformation functions
# Source: TransAnalyzerPipeline._sqrt, ._loga (trans_analyzer.py, lines 80-84)
# ============================================================================

def _sqrt(x):
    """Square root transformation: f(x) = sqrt(x)."""
    return np.sqrt(x)


def _loga(x, a=1):
    """Log transformation with pseudocount: f(x) = log(x + a)."""
    return np.log(x + a)


def _format_pseudocount(a_val):
    """
    Format a pseudocount value for use in layer names.

    Converts integer-valued floats to clean integers in the string
    representation. This ensures layer names like 'log(x+1)' instead
    of 'log(x+1.0)', which is critical for downstream matching against
    TRANSFORM_ORDER constants.

    Parameters
    ----------
    a_val : float
        The pseudocount value.

    Returns
    -------
    str
        Formatted string: '1' for 1.0, '10' for 10.0, '0.01' for 0.01.
    """
    if a_val == int(a_val):
        return str(int(a_val))
    return str(a_val)


def _acosh_transform(X, alpha=0.05):
    """
    Variance-stabilizing acosh transformation (Ahlmann-Eltze & Huber, 2023).

    Derives from the delta method applied to the NB mean-variance
    relationship Var[Y] = mu + alpha * mu^2:

        g(y) = (1 / sqrt(alpha)) * acosh(2 * alpha * y + 1)

    Applied to size-factor-normalized counts (y/s) where s is the simple
    size factor s_c = sum_g(y_gc) / L with L = mean of per-cell totals.

    Parameters
    ----------
    X : np.ndarray or scipy.sparse matrix
        Raw count matrix (cells x genes). Size factor normalization
        is applied internally.
    alpha : float
        NB overdispersion parameter (default 0.05, as recommended by
        the paper for a generic fixed value).

    Returns
    -------
    np.ndarray
        Transformed matrix (cells x genes), float32.
    """
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Compute simple size factors: s_c = total_c / mean(total)
    lib_sizes = X.sum(axis=1, keepdims=True)
    mean_lib = lib_sizes.mean()
    size_factors = lib_sizes / mean_lib

    # Size-factor-normalize
    Y = X / size_factors

    # Apply acosh VST: g(y) = (1/sqrt(alpha)) * acosh(2*alpha*y + 1)
    transformed = (1.0 / np.sqrt(alpha)) * np.arccosh(2.0 * alpha * Y + 1.0)

    return transformed.astype(np.float32)


def _pearson_residuals(X, theta=100):
    """
    Analytic Pearson residuals (Lause, Berens & Kobak, 2021).

    Computes residuals without per-gene regression, using a fixed
    overdispersion parameter theta:

        mu_ij = n_i * p_j
        r_ij  = (x_ij - mu_ij) / sqrt(mu_ij + mu_ij^2 / theta)

    where n_i is the library size of cell i and p_j = sum_i(x_ij) / sum_i(n_i)
    is the global gene frequency.

    Residuals are clipped to [-sqrt(n_cells), sqrt(n_cells)] to limit
    outlier influence, following Hafemeister & Satija (2019).

    Parameters
    ----------
    X : np.ndarray or scipy.sparse matrix
        Raw count matrix (cells x genes). Must be non-negative.
    theta : float
        Fixed NB overdispersion parameter (default 100, as used by
        Lause et al. 2021 and scTransform v2).

    Returns
    -------
    np.ndarray
        Clipped Pearson residuals (cells x genes), float32.
    """
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Library sizes and gene sums
    sums_cells = X.sum(axis=1, keepdims=True)   # (cells, 1)
    sums_genes = X.sum(axis=0, keepdims=True)   # (1, genes)
    sum_total = sums_genes.sum()

    # Expected counts under null model
    mu = sums_cells @ sums_genes / sum_total  # (cells, genes)

    # Pearson residuals with NB variance
    residuals = (X - mu) / np.sqrt(mu + mu ** 2 / theta)

    # Clip to [-sqrt(n), sqrt(n)]
    clip_val = np.sqrt(X.shape[0])
    np.clip(residuals, -clip_val, clip_val, out=residuals)

    return residuals.astype(np.float32)


def apply_standalone_transforms(raw_counts, alpha=0.05, theta=100):
    """
    Apply standalone transformations that operate on raw counts with
    internal size-factor handling, independent of any normalization.

    These are added to the benchmark alongside the per-normalization
    transformations (log, sqrt) but use a different computational
    approach: they compute their own size factors or expected counts
    internally rather than relying on an external normalization step.

    Parameters
    ----------
    raw_counts : np.ndarray or scipy.sparse matrix
        Raw count matrix (cells x genes).
    alpha : float
        NB overdispersion for acosh (default 0.05).
    theta : float
        NB overdispersion for Pearson residuals (default 100).
        Note: alpha and theta are related but NOT reciprocals here,
        because they follow different conventions in their respective
        papers. alpha is Var = mu + alpha*mu^2 (Ahlmann-Eltze),
        theta is Var = mu + mu^2/theta (Lause/scanpy).
        alpha=0.05 and theta=100 both represent mild overdispersion.

    Returns
    -------
    dict
        Mapping of {layer_name: transformed_matrix}.
    """
    results = {}

    print(f"  Transforming: none-acosh (alpha={alpha})")
    results['none-acosh'] = _acosh_transform(raw_counts, alpha=alpha)

    print(f"  Transforming: none-pearson_residuals (theta={theta})")
    results['none-pearson_residuals'] = _pearson_residuals(raw_counts, theta=theta)

    return results


# ============================================================================
# Phase 1: Per-layer transformation and statistics
# Source: TransAnalyzerPipeline.apply_transformations_single_layer (line 199)
#         TransAnalyzerPipeline.process_pca_and_stats_single_derived_layer (line 299)
# ============================================================================

def apply_transformations(base_data, base_layer, transformation_methods, pseudocounts):
    """
    Apply transformation functions to a single normalized expression matrix.

    For each transformation method, generates a derived matrix and stores it
    with a compound name (e.g. 'cpm-log(x+1)', 'deseq2-sqrt'). Methods like
    scTransform and SpaNorm are skipped because they handle transformation
    internally.

    Parameters
    ----------
    base_data : np.ndarray or pd.DataFrame
        The normalized expression matrix (cells x genes).
    base_layer : str
        Name of the base normalization layer (e.g. 'cpm', 'tmm').
    transformation_methods : list of str
        Transformations to apply. Recognized values:
        'none' (identity), 'sqrt', 'loga' (log with pseudocounts),
        'scTransform', 'spanorm' (both skipped).
    pseudocounts : list of float
        Pseudocount values for the log transformation (e.g. [0.01, 0.1, 0.5, 1, 10]).

    Returns
    -------
    dict
        Mapping of {derived_layer_name: transformed_matrix}.
        For 'none', the key is the base_layer name itself.
    """
    results = {}

    for trans in transformation_methods:
        if trans == 'none':
            # Identity: keep the original normalized data
            results[base_layer] = base_data

        elif trans == 'sqrt':
            key = f"{base_layer}-sqrt"
            print(f"  Transforming: {key}")
            results[key] = _sqrt(base_data)

        elif trans == 'loga':
            for a_val in pseudocounts:
                a_str = _format_pseudocount(a_val)
                key = f"{base_layer}-log(x+{a_str})"
                print(f"  Transforming: {key}")
                results[key] = _loga(base_data, a_val)

        elif trans in ('scTransform', 'spanorm'):
            # These methods handle transformation internally;
            # their outputs are already variance-stabilized
            continue

    return results


def compute_mean_variance(matrix):
    """
    Compute per-gene mean and variance from an expression matrix.

    Source: TransAnalyzerPipeline._calculate_mean_variance
            (trans_analyzer.py, line 86)

    Parameters
    ----------
    matrix : np.ndarray or scipy.sparse matrix
        Expression matrix (cells x genes).

    Returns
    -------
    tuple of (np.ndarray, np.ndarray)
        (mean_per_gene, variance_per_gene), each of shape (n_genes,).
    """
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    mean = np.asarray(matrix.mean(axis=0)).ravel()
    variance = np.asarray(matrix.var(axis=0)).ravel()
    return mean, variance


def compute_pca_and_stats(sd_obj, layer_name, data_matrix, dataset_name, technology):
    """
    Compute mean/variance and PCA for a single derived (transformed) matrix.

    Temporarily registers the matrix as a layer in the AnnData object,
    runs scanpy PCA, extracts loadings, and cleans up. Returns a DataFrame
    with per-gene statistics suitable for downstream aggregation.

    Source: TransAnalyzerPipeline.process_pca_and_stats_single_derived_layer
            (trans_analyzer.py, line 299)

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
        SpatialData object.
    layer_name : str
        Name for this derived layer (e.g. 'cpm-log(x+1)').
    data_matrix : np.ndarray or pd.DataFrame
        The transformed expression matrix (cells x genes).
    dataset_name : str
        Name of the dataset.
    technology : str
        iST platform name.

    Returns
    -------
    pd.DataFrame
        Per-gene statistics with columns: geneID, mean, variance,
        PC1, PC2, PC3, PC4, PC5, layer, dataset, technology.
        If PCA fails, the PC columns are absent.
    """
    adata = sd_obj.tables['filtered']
    gene_names = adata.var_names.to_numpy()

    # Temporarily register this matrix as a layer so scanpy can access it
    adata.layers[layer_name] = data_matrix

    # Compute per-gene mean and variance
    mean, variance = compute_mean_variance(data_matrix)

    df_stats = pd.DataFrame({
        "geneID": gene_names,
        "mean": mean,
        "variance": variance,
        "layer": layer_name,
        "dataset": dataset_name,
        "technology": technology,
    })

    # PCA: 5 components
    print(f"    Computing PCA for {layer_name}...")
    try:
        sc.pp.pca(adata, n_comps=5, layer=layer_name)

        pc_array = adata.varm['PCs']
        df_pcs = pd.DataFrame(
            pc_array,
            index=gene_names,
            columns=[f"PC{i+1}" for i in range(pc_array.shape[1])]
        ).reset_index().rename(columns={'index': 'geneID'})

        # Merge PCA loadings with the stats DataFrame
        df_stats = pd.merge(df_stats, df_pcs, on="geneID", how="inner")

        # Cleanup PCA artifacts from AnnData to keep it clean for the next layer
        del adata.uns['pca']
        del adata.varm['PCs']
        del adata.obsm['X_pca']

        print(f"    PCA complete: {pc_array.shape[1]} components")

    except Exception as e:
        print(f"    PCA failed for {layer_name}: {e}")

    # Remove the temporary layer
    if layer_name in adata.layers:
        del adata.layers[layer_name]

    return df_stats


# ============================================================================
# Phase 2: Aggregated variance stabilization metrics
# Source: 02_transformation.ipynb, cell 48 (lines 187-292)
# ============================================================================

def _compute_binned_cv(mean_vals, var_vals, n_bins=20):
    """
    Bin genes by mean expression into equal-count bins, compute median
    variance per bin, return CV across bins.

    Implements Equation 7 from the BISTRO manuscript:
        CV_var = std(median_var_per_bin) / mean(median_var_per_bin)

    A CV near zero indicates uniform variance across the expression range
    (successful stabilization), while high values reveal residual
    mean-variance dependence.

    Source: 02_transformation.ipynb, cell 48, line 187.

    Parameters
    ----------
    mean_vals : np.ndarray
        Per-gene mean expression values.
    var_vals : np.ndarray
        Per-gene variance values.
    n_bins : int, optional
        Number of equal-count bins (default 20).

    Returns
    -------
    float
        Coefficient of variation of median variance across bins.
        Returns np.nan if fewer than 3 valid bins.
    """
    order = np.argsort(mean_vals)
    var_sorted = var_vals[order]

    # Equal-count bins (each bin contains the same number of genes)
    bin_edges = np.array_split(var_sorted, n_bins)
    bin_medians = np.array([np.median(b) for b in bin_edges if len(b) > 0])

    # Remove any zero or negative medians (would distort CV)
    bin_medians = bin_medians[bin_medians > 0]
    if len(bin_medians) < 3:
        return np.nan

    return np.std(bin_medians) / np.mean(bin_medians)


def _compute_lowess_r2(mean_vals, var_vals, frac=0.6):
    """
    Fit LOWESS to mean vs variance, return R-squared of the smooth fit.

    R-squared near 0 indicates no mean-variance dependence (stabilized),
    while R-squared near 1 indicates strong residual dependence.

    Source: 02_transformation.ipynb, cell 48, line 207.

    Parameters
    ----------
    mean_vals : np.ndarray
        Per-gene mean expression values.
    var_vals : np.ndarray
        Per-gene variance values.
    frac : float, optional
        Smoothing fraction for LOWESS (default 0.6).

    Returns
    -------
    float
        R-squared of the LOWESS fit. Returns np.nan on failure.
    """
    try:
        fitted = sm.nonparametric.lowess(
            endog=var_vals, exog=mean_vals, frac=frac, return_sorted=True
        )
        # Interpolate predicted values at original x positions
        y_pred = np.interp(mean_vals, fitted[:, 0], fitted[:, 1])

        ss_res = np.sum((var_vals - y_pred) ** 2)
        ss_tot = np.sum((var_vals - np.mean(var_vals)) ** 2)

        if ss_tot == 0:
            return np.nan
        return 1.0 - ss_res / ss_tot
    except Exception:
        return np.nan


def compute_all_metrics(mean_var_stats, n_bins=20, lowess_frac=0.6):
    """
    Compute all four variance stabilization metrics per dataset x layer.

    For each (layer, dataset) group in the input DataFrame:
      1. Log-log slope: linear regression of log10(var) on log10(mean).
         Slope of 0 = stabilized, 1 = Poisson, 2 = Negative Binomial.
      2. Spearman rho: rank correlation between mean and variance.
         0 = stabilized.
      3. LOWESS R-squared: fraction of variance explained by a smooth
         mean-variance curve. 0 = stabilized.
      4. Binned CV (Equation 7): coefficient of variation of median
         variance across equal-count bins. 0 = stabilized.

    Layers belonging to RAW_ONLY_NORMS (scTransform, SpaNorm) are only
    evaluated in their 'raw' (untransformed) form, as they handle
    transformation internally.

    Source: 02_transformation.ipynb, cell 48, line 228.

    Parameters
    ----------
    mean_var_stats : pd.DataFrame
        Per-gene stats with columns: mean, variance, layer, dataset.
    n_bins : int, optional
        Number of equal-count bins for binned CV (default 20).
    lowess_frac : float, optional
        Smoothing fraction for LOWESS (default 0.6).

    Returns
    -------
    pd.DataFrame
        One row per (dataset, normalization, transformation) combination,
        with columns: normalization, transformation, dataset, slope,
        spearman, lowess_r2, binned_cv.
    """
    records = []

    for (layer, dataset), grp in mean_var_stats.groupby(["layer", "dataset"]):
        norm, trans = parse_layer(layer)

        # Skip spanorm base (its outputs are handled via compound names)
        if norm == "spanorm":
            continue
        # Skip transformed versions of methods that handle transformation internally
        if norm in RAW_ONLY_NORMS and trans != "raw":
            continue

        # Clean data: remove NaN, Inf, and non-positive values
        sub = grp[["mean", "variance"]].dropna()
        sub = sub[np.isfinite(sub["mean"]) & np.isfinite(sub["variance"])]
        sub = sub[(sub["mean"] > 0) & (sub["variance"] > 0)]

        if len(sub) < 20:
            print(f"  Skipping {layer} / {dataset}: only {len(sub)} valid genes")
            continue

        mean_vals = sub["mean"].values
        var_vals = sub["variance"].values

        # 1. Log-log slope: log10(var) ~ log10(mean)
        log_m, log_v = np.log10(mean_vals), np.log10(var_vals)
        mask = np.isfinite(log_m) & np.isfinite(log_v)
        try:
            slope, _, _, _, _ = linregress(log_m[mask], log_v[mask])
        except Exception:
            slope = np.nan

        # 2. Spearman rho(mean, variance)
        try:
            rho, _ = spearmanr(mean_vals, var_vals)
        except Exception:
            rho = np.nan

        # 3. LOWESS R-squared
        r2 = _compute_lowess_r2(mean_vals, var_vals, frac=lowess_frac)

        # 4. Binned CV (Equation 7)
        bcv = _compute_binned_cv(mean_vals, var_vals, n_bins=n_bins)

        records.append({
            "normalization":  norm,
            "transformation": trans,
            "dataset":        dataset,
            "slope":          slope,
            "spearman":       rho,
            "lowess_r2":      r2,
            "binned_cv":      bcv,
        })

    df = pd.DataFrame(records)
    print(f"Computed metrics for {len(df)} dataset x layer combinations")
    print(f"  Datasets:        {df['dataset'].nunique()}")
    print(f"  Normalizations:  {sorted(df['normalization'].unique())}")
    print(f"  Transformations: {sorted(df['transformation'].unique())}")
    return df


# ============================================================================
# Phase 2b: PC1-mean correlation
# Source: 02_transformation.ipynb, cell 46 (line 122)
# ============================================================================

def compute_pc1_correlations(mean_var_stats, pc_col="PC1"):
    """
    Compute Spearman r(mean, |PC1|) per dataset x layer.

    A high absolute correlation indicates that PC1 is dominated by
    expression magnitude rather than capturing biological variance,
    reflecting a failure of the transformation to decouple the first
    principal component from the expression level.

    Implements the analysis described in Section 2.3 and referenced
    in Equation 6 of the BISTRO manuscript.

    Source: 02_transformation.ipynb, cell 46, line 122.

    Parameters
    ----------
    mean_var_stats : pd.DataFrame
        Per-gene stats with columns: mean, PC1, layer, dataset.
    pc_col : str, optional
        Column name for the PC1 loading (default 'PC1').

    Returns
    -------
    pd.DataFrame
        One row per (dataset, normalization, transformation) with
        columns: normalization, transformation, dataset, correlation.
    """
    records = []

    for (layer, dataset), grp in mean_var_stats.groupby(["layer", "dataset"]):
        sub = grp[["mean", pc_col]].dropna()
        sub = sub[np.isfinite(sub["mean"]) & np.isfinite(sub[pc_col])]
        if len(sub) < 10:
            continue

        norm, trans = parse_layer(layer)

        # Skip compound SpaNorm with transformations (already handled internally)
        if norm == "spanorm" or ("spanorm" in norm and trans != "raw"):
            continue
        # Skip transformed scTransform (already variance-stabilized)
        if norm == "scTransform" and trans != "raw":
            continue

        try:
            r, _ = spearmanr(sub["mean"].values, np.abs(sub[pc_col].values))
            records.append({
                "normalization": norm,
                "transformation": trans,
                "dataset": dataset,
                "correlation": r,
            })
        except Exception:
            continue

    df = pd.DataFrame(records)
    if not df.empty:
        print(f"Computed PC1 correlations for {len(df)} combinations")
    return df


# ============================================================================
# Pipeline: run transformation analysis end-to-end
# Source: 02_run_trans_analyzer.py main loop (lines 248-293)
# ============================================================================

def run_transformation_pipeline(
    zarr_path, nextflow_output, output_dir, dataset_name=None,
    technology=None, pseudocounts=None, alpha=0.05, theta=100,
):
    """
    Run the full transformation analysis across all normalization layers.

    This function orchestrates the complete analysis:
    1. Load the SpatialData object
    2. Discover normalization layers from Nextflow output
    3. Phase 1: For each normalization layer:
       a. Load the normalized expression matrix
       b. Apply all transformations (sqrt, log with pseudocounts, identity)
       c. Compute per-gene mean, variance, and PCA (5 components)
       d. Unload to free memory
    4. Standalone transforms on raw counts:
       a. Acosh VST (Ahlmann-Eltze & Huber, 2023)
       b. Analytic Pearson residuals (Lause, Berens & Kobak, 2021)
    5. Phase 2: Aggregate metrics across all combinations
       a. Compute log-log slope, Spearman rho, LOWESS R-squared, binned CV
       b. Compute PC1-mean correlations
    6. Save all outputs

    Parameters
    ----------
    zarr_path : str
        Path to the QC'd .zarr file.
    nextflow_output : str
        Path to the Nextflow output directory (containing norm/ subdirectory).
    output_dir : str
        Directory to write output files.
    dataset_name : str, optional
        Name for the dataset. If None, inferred from the zarr filename.
    technology : str, optional
        iST platform name ('CosMx', 'Xenium', 'MERFISH'). If None, set to 'unknown'.
    pseudocounts : list of float, optional
        Pseudocount values for log transformation (default [0.01, 0.1, 0.5, 1, 10]).
    alpha : float, optional
        NB overdispersion for acosh transform (default 0.05, per
        Ahlmann-Eltze & Huber 2023). Parameterization: Var = mu + alpha*mu^2.
    theta : float, optional
        NB overdispersion for Pearson residuals (default 100, per
        Lause et al. 2021). Parameterization: Var = mu + mu^2/theta.
    """
    os.makedirs(output_dir, exist_ok=True)

    if dataset_name is None:
        dataset_name = os.path.basename(zarr_path).replace('.zarr', '')
    if technology is None:
        technology = "unknown"
    if pseudocounts is None:
        pseudocounts = [0.01, 0.1, 0.5, 1, 10]

    # Transformation methods to apply to each scale-factor-based normalization
    transforms = ['sqrt', 'loga', 'scTransform', 'none', 'spanorm']

    # ---- Load data ----
    print(f"\n{'='*60}")
    print(f"BISTRO Transformation Analysis: {dataset_name}")
    print(f"{'='*60}")

    sd_obj = load_spatialdata(zarr_path)
    norm_map = discover_norm_layers(nextflow_output)

    if not norm_map:
        print("No normalization layers found. Exiting.")
        return

    # ---- Phase 1: Per-layer statistics ----
    print(f"\n{'='*60}")
    print("Phase 1: Per-layer transformation and statistics")
    print(f"{'='*60}")

    all_stats_dfs = []

    for layer in sorted(norm_map.keys()):
        # Skip spanorm base layer (it is handled via compound names)
        if layer == 'spanorm':
            continue

        print(f"\n--- Processing base layer: {layer} ---")
        try:
            # Load normalized expression matrix
            load_layer_into_spatialdata(sd_obj, layer, norm_map)
            base_data = sd_obj.tables['filtered'].layers[layer]

            # Apply all transformations to this base layer
            derived_layers = apply_transformations(
                base_data, layer, transforms, pseudocounts
            )

            # Compute mean/variance/PCA for each derived layer
            for derived_name, data_matrix in derived_layers.items():
                print(f"  Analyzing: {derived_name}")
                df_stats = compute_pca_and_stats(
                    sd_obj, derived_name, data_matrix, dataset_name, technology
                )
                all_stats_dfs.append(df_stats)
                del data_matrix

            # Clean up
            del derived_layers
            unload_layer(sd_obj, layer)
            gc.collect()

        except Exception as e:
            print(f"  FAILED: {e}")
            unload_layer(sd_obj, layer)
            continue

    # ---- Standalone transforms (applied once to raw counts) ----
    # These transformations compute their own size factors / expected counts
    # internally and do not depend on any external normalization step.
    # They appear as 'none-acosh' and 'none-pearson_residuals' in the output.
    print(f"\n{'='*60}")
    print("Standalone transforms on raw counts")
    print(f"{'='*60}")
    try:
        raw_X = sd_obj.tables['filtered'].X
        standalone_layers = apply_standalone_transforms(
            raw_X, alpha=alpha, theta=theta
        )
        for derived_name, data_matrix in standalone_layers.items():
            print(f"  Analyzing: {derived_name}")
            df_stats = compute_pca_and_stats(
                sd_obj, derived_name, data_matrix, dataset_name, technology
            )
            all_stats_dfs.append(df_stats)
            del data_matrix
        del standalone_layers
        gc.collect()
    except Exception as e:
        print(f"  FAILED standalone transforms: {e}")

    if not all_stats_dfs:
        print("No results generated.")
        return

    # Concatenate all per-gene stats into a single DataFrame
    mean_var_stats = pd.concat(all_stats_dfs, ignore_index=True)

    # Save Phase 1 output
    stats_path = os.path.join(output_dir, f"{dataset_name}_mean_variance_stats.csv")
    mean_var_stats.to_csv(stats_path, index=False)
    print(f"\nSaved per-gene stats: {stats_path}")
    print(f"  Shape: {mean_var_stats.shape}")
    print(f"  Layers: {mean_var_stats['layer'].nunique()}")

    # ---- Phase 2: Aggregated metrics ----
    print(f"\n{'='*60}")
    print("Phase 2: Aggregated variance stabilization metrics")
    print(f"{'='*60}")

    df_metrics = compute_all_metrics(mean_var_stats)
    metrics_path = os.path.join(output_dir, f"{dataset_name}_transformation_summary.csv")
    df_metrics.to_csv(metrics_path, index=False)
    print(f"Saved transformation summary: {metrics_path}")

    # ---- Phase 2b: PC1 correlation ----
    print(f"\n{'='*60}")
    print("Phase 2b: PC1-mean expression correlation")
    print(f"{'='*60}")

    if 'PC1' in mean_var_stats.columns:
        df_corr = compute_pc1_correlations(mean_var_stats)
        corr_path = os.path.join(output_dir, f"{dataset_name}_pc1_correlation.csv")
        df_corr.to_csv(corr_path, index=False)
        print(f"Saved PC1 correlations: {corr_path}")
    else:
        print("Warning: PC1 column not found in mean_var_stats, "
              "skipping correlation analysis.")

    print(f"\nTransformation analysis complete for {dataset_name}.")


# ============================================================================
# CLI
# ============================================================================

def build_parser():
    """Build argument parser for the transformation analysis CLI."""
    p = argparse.ArgumentParser(
        description="BISTRO: Transformation / Variance Stabilization Analysis. "
                    "Applies transformations to each normalization output and "
                    "computes mean-variance metrics, PCA loadings, and "
                    "variance stabilization summaries.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--zarr", required=True,
                   help="Path to the QC'd .zarr file")
    p.add_argument("--nextflow_output", required=True,
                   help="Path to the Nextflow output directory "
                        "(must contain a norm/ subdirectory)")
    p.add_argument("--output_dir", required=True,
                   help="Directory to write output files")
    p.add_argument("--dataset_name", default=None,
                   help="Dataset name (default: inferred from zarr path)")
    p.add_argument("--technology", default=None,
                   help="iST platform (CosMx, Xenium, MERFISH)")
    p.add_argument("--pseudocounts", type=float, nargs="+",
                   default=[0.01, 0.1, 0.5, 1, 10],
                   help="Pseudocount values for log transformation")
    p.add_argument("--alpha", type=float, default=0.05,
                   help="NB overdispersion for acosh VST "
                        "(default 0.05, per Ahlmann-Eltze & Huber 2023)")
    p.add_argument("--theta", type=float, default=100,
                   help="NB overdispersion for analytic Pearson residuals "
                        "(default 100, per Lause et al. 2021)")
    return p


def main():
    """CLI entry point."""
    args = build_parser().parse_args()
    run_transformation_pipeline(
        zarr_path=args.zarr,
        nextflow_output=args.nextflow_output,
        output_dir=args.output_dir,
        dataset_name=args.dataset_name,
        technology=args.technology,
        pseudocounts=args.pseudocounts,
        alpha=args.alpha,
        theta=args.theta,
    )


if __name__ == "__main__":
    main()