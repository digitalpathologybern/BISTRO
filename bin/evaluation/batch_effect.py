#!/usr/bin/env python
"""
BISTRO Batch Effect Evaluation
==============================

Evaluates the intra-slide FOV-based batch effect for each normalization
method by fitting an OLS model (tissue only) and a Mixed-Effects Linear
Model (MELM) with FOV as a random intercept.

Implements Equations 1-5 from the BISTRO manuscript (Section 4.6).

Outputs
-------
  {dataset}_batch_effect_summary.csv
      One row per normalization method with: AIC, BIC for both models,
      delta-BIC, LRT p-value, random intercept variance with CI,
      drift slope, drift Pearson r and p-value.

  {dataset}_random_intercepts.csv
      Per-FOV random intercepts across all normalizations, including
      the relative bias b_j = exp(u_j) - 1 (Equation 4).

  {dataset}_batch_effect_models.pkl
      Serialized full statsmodels model objects (OLS, MELM) for every
      normalization method, enabling downstream inspection or re-analysis.

  {dataset}_fov_summary.csv
      Per-FOV mean library size with spatial coordinates (fov_center_x_um,
      fov_center_y_um), computed on raw counts. Used for the rasterized
      FOV heatmap in the HTML report.

  {dataset}_tissue_ls.csv
      Per-tissue-region mean and standard deviation of library size,
      computed on raw counts. Used for the tissue barplot in the HTML report.

Usage
-----
    python batch_effect.py \\
        --zarr /path/to/dataset.zarr \\
        --nextflow_output /path/to/nextflow_output/ \\
        --technology CosMx \\
        --tissue_annotation /path/to/annotations.csv \\
        --output_dir /path/to/output/
"""

import argparse
import gc
import os
import pickle
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import pearsonr

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'utils'))
from helpers import calculate_library_size, var_ci_chisq, lrt_pvalue, bootstrap_var_ci
from data_loader import (
    load_spatialdata, discover_norm_layers, load_layer_into_spatialdata,
    unload_layer, load_fov_from_metadata, load_tissue_annotations,
)


# ============================================================================
# Core analysis: process a single normalization layer
# Source: NormAnalyzerPipeline._process_layer_fast
#         (norm_analyzer.py, lines 1631-1734)
# ============================================================================

def evaluate_batch_effect_single_layer(
    sd_obj,
    layer_name,
    dataset_name,
    technology,
    libsize_column='library_size',
    tissue_column='tissue_annotations',
    fov_column='fov',
    use_log=True,
    n_jobs=None,
    ci_scheme="cluster",
    n_boot_ci=200,
    vc_column=None,
):
    """
    Fit OLS and MELM models for a single normalization layer to evaluate
    the FOV-based batch effect.

    This function implements the statistical framework described in
    Section 4.6 (Batch effect evaluation) of the BISTRO manuscript:

      - Equation 1 (LR model):
            y_i = b0 + sum_k(bk * I[tissue_i = k]) + e_i
        Fitted via OLS with tissue region as the sole fixed effect.

      - Equation 2 (MELM):
            y_i = b0 + sum_k(bk * I[tissue_i = k]) + u_j[i] + e_i
        where u_j ~ N(0, tau^2) is the FOV-specific random intercept.
        Fitted twice: once with ML (for BIC/AIC/LRT model comparison)
        and once with REML (for unbiased random intercept extraction
        and drift analysis).

      - Equation 3 (Likelihood ratio test):
            Lambda = 2 * (loglik_MELM - loglik_LR)
        tested against the boundary mixture 0.5*chi2_1 + 0.5*chi2_0, since
        the FOV variance cannot be negative. Uses the ML fits (consistent
        with BIC comparison).

      - With vc_column (the TMA's patient), FOVs are nested in it:
            y_i = b0 + sum_k(bk * I[tissue_i = k]) + v_p[i] + u_j[i] + e_i
        v_p ~ N(0, sigma_p^2) is the patient intercept shared by that
        patient's FOVs and u_j ~ N(0, tau^2) the FOV intercept within it.
        The LRT then tests tau^2 against the model that keeps v_p.

      - Equation 4 (Relative bias per FOV):
            b_j = exp(u_hat_j) - 1
        Extracted from the REML fit (unbiased variance components).

      - Equation 5 (Bootstrap CI for tau^2):
        Cluster bootstrap: resample the top-level groups (FOVs, or patients
        with their FOVs when nested), refit REML MELM, collect tau-squared,
        take quantiles. Replaces the chi-squared approximation which has
        poor coverage for small group counts.

    Additionally, this function computes a drift regression by fitting
    the random intercepts against the FOV acquisition index using OLS,
    to quantify systematic signal degradation during imaging.

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
        SpatialData object with 'filtered' table. The specified layer must
        already be loaded into sd_obj.tables['filtered'].layers[layer_name].
    layer_name : str
        Name of the normalization layer to evaluate (e.g. 'cpm', 'scTransform').
    dataset_name : str
        Name of the dataset, used for labeling output rows.
    technology : str
        iST platform name ('CosMx', 'Xenium', 'MERSCOPE').
    libsize_column : str, optional
        Column name for library size in obs (default 'library_size').
    tissue_column : str, optional
        Column name for tissue region annotations in obs
        (default 'tissue_annotations').
    fov_column : str, optional
        Column name for FOV identifier in obs (default 'fov').
    use_log : bool, optional
        If True (default), use log1p(library_size) as the dependent variable
        in the regression models.
    n_boot_ci : int, optional
        Number of bootstrap iterations for the tau-squared CI (default 200).

    Returns
    -------
    dict
        Dictionary containing model comparison statistics (from ML fit),
        random intercepts and drift analysis (from REML fit), and
        bootstrap CI for the random intercept variance.
    """
    print(f"Evaluating batch effect for layer: {layer_name}")

    matrix = sd_obj.tables['filtered'].layers[layer_name]
    adata = sd_obj.tables['filtered']

    # Compute library size from the normalized expression matrix
    library_sizes = calculate_library_size(matrix)
    adata.obs['library_size'] = library_sizes
    adata.obs['log_LS'] = np.log1p(library_sizes)

    n_nan = int(np.isnan(adata.obs['log_LS']).sum())
    n_neg = int((library_sizes < 0).sum())

    # True only when the row sum is a non-negative per-cell total.
    response_is_library_size = (n_nan == 0 and n_neg == 0)

    if not response_is_library_size and use_log:
        print(f"  NOTE: layer '{layer_name}' has {n_neg} cells with a negative "
              f"row sum and {n_nan} non-finite log values. Its row sum is a "
              f"residual total, not a library size; falling back to the "
              f"untransformed response and marking the library-size drift and "
              f"relative-bias statistics undefined for this layer.")
        use_log = False

    # Prepare modelling dataframe with required columns
    nested = bool(vc_column and vc_column in adata.obs.columns)
    required_cols = [libsize_column, tissue_column, fov_column]
    if use_log:
        required_cols.append('log_LS')
    if nested:
        required_cols.append(vc_column)
    df = adata.obs[required_cols].copy()
    # A blank label is a missing label, not a level of its own.
    for col in [tissue_column] + ([vc_column] if nested else []):
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].astype(object).where(
                df[col].astype(str).str.strip() != '', np.nan)
    df.dropna(inplace=True)

    df[tissue_column] = df[tissue_column].astype('category')
    df[fov_column] = df[fov_column].astype('int')
    df['fov_complete'] = 'FOV' + df[fov_column].astype(str)
    fov_col = 'fov_complete'

    # Without an extra variance component FOV is the grouping factor. With one
    # (the TMA's patient) FOVs are nested in it: the patient is the group and the
    # FOV a variance component inside it, so the cores of one patient share that
    # patient's intercept and the FOV variance is what remains within a patient.
    if nested:
        df[vc_column] = df[vc_column].astype(str).astype('category')
        n_per_fov = df.groupby(fov_col)[vc_column].nunique()
        if (n_per_fov > 1).any():
            raise ValueError(
                f"{int((n_per_fov > 1).sum())} FOVs span more than one "
                f"'{vc_column}'; FOVs must be nested in '{vc_column}'.")
        group_col = vc_column
        vc_formula = {'fov': f"0 + C({fov_col})"}
        print(f"  Nested design: {df[fov_col].nunique()} FOVs in "
              f"{df[vc_column].nunique()} levels of '{vc_column}'")
    else:
        group_col = fov_col
        vc_formula = None

    y_col = 'log_LS' if use_log else libsize_column
    formula = f"{y_col} ~ C({tissue_column})"

    # ---- Model 1: OLS with tissue only (Equation 1) ----
    print(f"  Fitting OLS: {formula}")
    model_ols = smf.ols(formula=formula, data=df).fit()

    # ---- Model 2: OLS with tissue + FOV as fixed effect (diagnostic) ----
    formula_fov = f"{y_col} ~ C({tissue_column}) + C({fov_col})"
    print(f"  Fitting OLS: {formula_fov}")
    model_ols_fov = smf.ols(formula=formula_fov, data=df).fit()

    # ---- Model 3a: MELM with ML (for BIC/AIC/LRT model comparison) ----
    def _fit_melm(reml, vc=vc_formula):
        return smf.mixedlm(
            formula=formula, groups=df[group_col], data=df, re_formula='~1',
            vc_formula=vc,
        ).fit(method=["powell"], reml=reml)

    model_desc = (f"(1|{vc_column}) + (1|{vc_column}:{fov_col})" if nested
                  else f"(1|{fov_col})")
    print(f"  Fitting MELM (ML): {formula} + {model_desc}")
    model_melm_ml = _fit_melm(reml=False)

    # ---- Model 3b: MELM with REML (for unbiased intercept extraction) ----
    print(f"  Fitting MELM (REML): {formula} + {model_desc}")
    model_melm_reml = _fit_melm(reml=True)

    # The model without the FOV term, against which the FOV term is tested:
    # the tissue-only OLS, or with nesting the model that keeps the patient.
    if nested:
        print(f"  Fitting MELM (ML) without FOV: {formula} + (1|{vc_column})")
        model_no_fov = _fit_melm(reml=False, vc=None)
    else:
        model_no_fov = model_ols

    # ---- Compute AIC and BIC from ML fit (for model comparison) ----
    aic_ols, bic_ols = model_no_fov.aic, model_no_fov.bic
    aic_ols_fov, bic_ols_fov = model_ols_fov.aic, model_ols_fov.bic
    aic_melm, bic_melm = model_melm_ml.aic, model_melm_ml.bic

    print(f"  Without FOV:      AIC = {aic_ols:.1f}, BIC = {bic_ols:.1f}")
    print(f"  OLS (tissue+FOV): AIC = {aic_ols_fov:.1f}, BIC = {bic_ols_fov:.1f}")
    print(f"  MELM ML:          AIC = {aic_melm:.1f}, BIC = {bic_melm:.1f}")

    # ---- LRT p-value: without FOV vs MELM-ML (Equation 3) ----
    # One variance component, the FOV's, is tested in both designs.
    p_value = lrt_pvalue(model_no_fov, model_melm_ml, n_vc=1)
    print(f"  LRT p-value (FOV variance): {p_value:.4e}")

    # ---- Extract random intercepts from REML fit (Equation 4) ----
    # Using REML for intercepts because ML underestimates variance components,
    # which biases the per-FOV random intercepts used in visualization and drift.
    re_df = _fov_random_effects(model_melm_reml, vc_column if nested else None)
    if use_log and response_is_library_size:
        re_df['relative_bias'] = np.exp(re_df['random_intercept']) - 1
    else:
        re_df['relative_bias'] = np.nan
    re_df['layer'] = layer_name
    re_df['dataset'] = dataset_name

    # ---- Variance estimates from both fits ----
    # tau^2 is always the FOV variance; with nesting it is the variance of the
    # FOV component and the group variance is the patient's.
    var_vc_reml = np.nan
    if nested:
        var_hat_ml = model_melm_ml.vcomp[0]
        var_hat_reml = model_melm_reml.vcomp[0]
        var_vc_reml = model_melm_reml.cov_re.iloc[0, 0]
    else:
        var_hat_ml = model_melm_ml.cov_re.iloc[0, 0]
        var_hat_reml = model_melm_reml.cov_re.iloc[0, 0]

    print(f"  Var(tau_fov^2) ML   = {var_hat_ml:.6f}")
    print(f"  Var(tau_fov^2) REML = {var_hat_reml:.6f}")
    if nested:
        print(f"  Var({vc_column}) REML = {var_vc_reml:.6f}")

    # ---- Bootstrap CI for tau-squared (replaces chi-squared, Equation 5) ----
    print(f"  Computing bootstrap CI for tau^2 ({n_boot_ci} iterations)...")
    var_ci = bootstrap_var_ci(
        df, formula, group_col, n_boot=n_boot_ci, alpha=0.05, seed=42,
        n_jobs=n_jobs, scheme=ci_scheme, vc_formula=vc_formula,
        nested_col=fov_col if nested else None,
    )
    print(f"  Bootstrap CI = [{var_ci[0]:.6f}, {var_ci[1]:.6f}]")

    # ---- Drift regression: REML intercepts vs FOV acquisition index ----
    # Drift response: exp(u) - 1 for library-size layers, u otherwise.
    # Both the rank and the FOV-identifier axis are reported.
    from scipy.stats import linregress

    re_sorted = re_df.sort_values('fov')
    if use_log and response_is_library_size:
        intercepts = re_sorted['relative_bias'].to_numpy(dtype=float)
        drift_response = 'relative_bias'
    else:
        intercepts = re_sorted['random_intercept'].to_numpy(dtype=float)
        drift_response = 'random_intercept'
    axis_rank = np.arange(len(re_sorted), dtype=float)
    axis_id = re_sorted['fov'].to_numpy(dtype=float)

    finite = np.isfinite(intercepts)
    n_nonfinite = int((~finite).sum())
    if n_nonfinite:
        print(f"  WARNING: {n_nonfinite} non-finite random intercepts excluded "
              f"from the drift regression")

    def _drift(x):
        m = finite & np.isfinite(x)
        if m.sum() <= 2 or np.std(intercepts[m]) == 0 or np.std(x[m]) == 0:
            return (np.nan,) * 4
        res = linregress(x[m], intercepts[m])
        return res.rvalue, res.pvalue, res.slope, res.stderr

    drift_r, drift_p, drift_slope, drift_se = _drift(axis_rank)
    drift_r_id, drift_p_id, drift_slope_id, drift_se_id = _drift(axis_id)

    if not response_is_library_size:
        print(f"  Drift on '{layer_name}': computed on the residual response; "
              f"NOT a library-size drift and reported as undefined.")

    print(f"  Drift (rank axis, on {drift_response}): r = {drift_r:.4f}, "
          f"p = {drift_p:.4e}, slope = {drift_slope:.6g} +/- {drift_se:.3g}")

    return {
        'layer': layer_name,
        'dataset': dataset_name,
        'technology': technology,
        # Model comparison (from ML fit)
        'aic_no_fov': aic_ols,
        'bic_no_fov': bic_ols,
        'aic_with_fov': aic_ols_fov,
        'bic_with_fov': bic_ols_fov,
        'aic_mixedlm_random_intercept': aic_melm,
        'bic_mixedlm_random_intercept': bic_melm,
        'pvalue_ols_no_fov_vs_mixedlm_random_intercept': p_value,
        # Model objects
        'model_no_fov': model_no_fov,
        'model_with_fov': model_ols_fov,
        'mixedlm_model_ml': model_melm_ml,
        'mixedlm_model_reml': model_melm_reml,
        # Variance estimates
        'var_mixedlm_ml': var_hat_ml,
        'var_mixedlm_reml': var_hat_reml,
        'var_ci_lower': var_ci[0],
        'var_ci_upper': var_ci[1],
        'ci_method': f'bootstrap_{ci_scheme}',
        'n_boot_ci': n_boot_ci,
        'vc_column': vc_column if vc_column else '',
        'var_vc_reml': var_vc_reml,
        # Random intercepts (from REML fit)
        'random_effects_mixedlm_model_random_intercept': re_df,
        # Drift analysis (from REML intercepts)
        'drift_pearson_r': drift_r,
        'drift_pearson_p': drift_p,
        'drift_slope': drift_slope,
        'drift_slope_se': drift_se,
        'drift_axis': 'rank',
        'drift_response': drift_response,
        'drift_pearson_r_fovid': drift_r_id,
        'drift_pearson_p_fovid': drift_p_id,
        'drift_slope_fovid': drift_slope_id,
        'drift_slope_se_fovid': drift_se_id,
        'n_nonfinite_intercepts': n_nonfinite,
        # Interpretability flags. See the response classification above.
        'response_is_library_size': response_is_library_size,
        'response_scale': 'log1p' if use_log else 'linear',
        'response_sd': float(np.nanstd(df[y_col].to_numpy(dtype=float))),
        # Relative test: the response has no variance to model.
        'response_degenerate': bool(
            np.nanstd(df[y_col].to_numpy(dtype=float))
            / max(abs(float(np.nanmean(df[y_col].to_numpy(dtype=float)))), 1e-12)
            < 1e-6),
        'lrt_reference': '0.5*chi2_1+0.5*chi2_0',
        'model_structure': f'fov_in_{vc_column}' if nested else 'fov',
        'fixed_effect_column': tissue_column,
    }


def _fov_random_effects(model, vc_column=None):
    """
    Per-FOV random intercepts from a fitted MixedLM.

    Without nesting each group is a FOV. With nesting each group is a level of
    vc_column; its FOVs are the variance-component entries, named
    'fov[C(fov_complete)[FOV12]]', and the group's own intercept is carried
    alongside each of its FOVs.
    """
    rows = []
    for group, eff in model.random_effects.items():
        if vc_column is None:
            rows.append({'fov': group, 'random_intercept': eff['Group']})
            continue
        for name, value in eff.items():
            if name == 'Group':
                continue
            rows.append({
                'fov': name.rsplit('[', 1)[1].rstrip(']'),
                'random_intercept': value,
                vc_column: group,
                f'{vc_column}_intercept': eff['Group'],
            })
    re_df = pd.DataFrame(rows)
    re_df['fov'] = re_df['fov'].str.replace('FOV', '').astype(int)
    return re_df


# ============================================================================
# FOV and tissue summaries for reporting
# ============================================================================

def compute_fov_and_tissue_summaries(sd_obj, dataset_name,
                                     tissue_column='tissue_annotations',
                                     fov_column='fov'):
    """
    Compute per-FOV mean library size (with spatial coordinates) and
    per-tissue-region mean library size from the raw count matrix.

    These summaries are computed once on the unnormalized data and
    used to generate the rasterized FOV heatmap and the tissue region
    barplot in the HTML report.

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
        SpatialData object with raw counts in .X of the 'filtered' table.
    dataset_name : str
        Dataset name for labeling output rows.
    tissue_column : str, optional
        Column name for tissue region annotations (default 'tissue_annotations').
    fov_column : str, optional
        Column name for FOV identifier (default 'fov').

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        fov_summary : DataFrame with columns [fov, mean_ls, n_cells,
                      fov_center_x_um, fov_center_y_um, dataset].
        tissue_ls   : DataFrame with columns [tissue_region, mean_ls,
                      std_ls, n_cells, dataset]. Empty if tissue_column
                      is not present in obs.
    """
    adata = sd_obj.tables['filtered']

    # Compute raw library sizes from the unprocessed count matrix
    raw_ls = calculate_library_size(adata.X)
    adata.obs['raw_library_size'] = raw_ls

    # ---- Per-FOV mean library size with spatial coordinates ----
    fov_stats = (
        adata.obs
        .groupby(fov_column, observed=True)
        .agg(
            mean_ls=('raw_library_size', 'mean'),
            n_cells=('raw_library_size', 'count'),
        )
        .reset_index()
    )

    # Add FOV center coordinates if available (present for CosMx natively,
    # added by rasterize_fov_grid for Xenium/MERSCOPE)
    for coord_col in ['fov_center_x_um', 'fov_center_y_um']:
        if coord_col in adata.obs.columns:
            fov_coords = (
                adata.obs[[fov_column, coord_col]]
                .drop_duplicates()
            )
            fov_stats = fov_stats.merge(fov_coords, on=fov_column, how='left')

    fov_stats['dataset'] = dataset_name

    print(f"  FOV summary: {len(fov_stats)} FOVs, "
          f"mean LS range = [{fov_stats['mean_ls'].min():.1f}, {fov_stats['mean_ls'].max():.1f}]")

    # ---- Per-tissue-region mean library size ----
    tissue_ls = pd.DataFrame()
    if tissue_column in adata.obs.columns:
        tissue_ls = (
            adata.obs
            .groupby(tissue_column, observed=True)
            .agg(
                mean_ls=('raw_library_size', 'mean'),
                std_ls=('raw_library_size', 'std'),
                n_cells=('raw_library_size', 'count'),
            )
            .reset_index()
            .rename(columns={tissue_column: 'tissue_region'})
            .sort_values('mean_ls', ascending=False)
        )
        tissue_ls['dataset'] = dataset_name

        print(f"  Tissue summary: {len(tissue_ls)} regions")
        for _, row in tissue_ls.iterrows():
            print(f"    {row['tissue_region']}: mean={row['mean_ls']:.1f}, "
                  f"std={row['std_ls']:.1f}, n={row['n_cells']}")
    else:
        print(f"  Warning: tissue column '{tissue_column}' not found in obs, "
              f"skipping tissue summary.")

    return fov_stats, tissue_ls


# ============================================================================
# Perpendicular FOV drift comparison
# ============================================================================

def compute_drift_comparison(sd_obj, dataset_name, technology,
                             tissue_column='tissue_annotations',
                             use_log=True):
    """
    Compute the FOV signal drift for both the default (row-primary) and
    perpendicular (col-primary) FOV orderings on the raw count data.

    This is only meaningful for rasterized platforms (Xenium, MERSCOPE)
    where the 'fov_perp' column is present in the metadata. For CosMx
    (native FOVs), returns None.

    The function fits a MELM with random intercept for each FOV ordering,
    extracts the per-FOV random intercepts, and computes the drift
    regression (relative bias vs FOV acquisition index).

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
        SpatialData object with raw counts and FOV columns in obs.
    dataset_name : str
        Dataset name for labeling.
    technology : str
        iST platform name.
    tissue_column : str, optional
        Column name for tissue region annotations.
    use_log : bool, optional
        Whether to log-transform library sizes.

    Returns
    -------
    dict or None
        Dictionary with keys 'row' and 'col', each containing:
          - 'random_intercepts': pd.DataFrame with per-FOV intercepts
          - 'drift_r', 'drift_p', 'drift_slope': drift statistics
          - 'fov_centers': pd.DataFrame with FOV spatial coordinates
        Returns None if 'fov_perp' is not in obs.
    """
    adata = sd_obj.tables['filtered']

    if 'fov_perp' not in adata.obs.columns:
        print("  No perpendicular FOV column found (native FOVs). Skipping drift comparison.")
        return None

    print(f"\n{'='*60}")
    print("Computing drift comparison: row-primary vs col-primary FOV ordering")
    print(f"{'='*60}")

    raw_ls = calculate_library_size(adata.X)
    adata.obs['library_size'] = raw_ls
    adata.obs['log_LS'] = np.log1p(raw_ls)

    y_col = 'log_LS' if use_log else 'library_size'
    results = {}

    for fov_col, label in [('fov', 'row'), ('fov_perp', 'col')]:
        print(f"\n  --- {label}-primary scan (column: {fov_col}) ---")

        df = adata.obs[[y_col, tissue_column, fov_col]].copy()
        df.dropna(inplace=True)
        df[tissue_column] = df[tissue_column].astype('category')
        df[fov_col] = df[fov_col].astype('int')
        df['fov_str'] = 'FOV' + df[fov_col].astype(str)
        formula = f"{y_col} ~ C({tissue_column})"

        # Fit ML for BIC comparison
        model_ml = smf.mixedlm(
            formula=formula, groups=df['fov_str'], data=df, re_formula='~1'
        ).fit(method=["powell"], reml=False)

        # Fit REML for unbiased intercept extraction
        model_reml = smf.mixedlm(
            formula=formula, groups=df['fov_str'], data=df, re_formula='~1'
        ).fit(method=["powell"], reml=True)

        # Extract random intercepts from REML fit
        re_df = pd.DataFrame(model_reml.random_effects).T.reset_index()
        re_df = re_df.rename(columns={'index': 'fov', 'Group': 'random_intercept'})
        re_df['fov'] = re_df['fov'].str.replace('FOV', '').astype(int)
        re_df['relative_bias'] = np.exp(re_df['random_intercept']) - 1
        re_df['scan_direction'] = label
        re_df['dataset'] = dataset_name

        # Drift regression on REML intercepts
        re_sorted = re_df.sort_values('fov')
        fov_idx = np.arange(len(re_sorted))
        bias = re_sorted['relative_bias'].values

        if len(fov_idx) > 2 and np.std(bias) > 0:
            drift_r, drift_p = pearsonr(fov_idx, bias)
            from sklearn.linear_model import LinearRegression
            lr = LinearRegression().fit(fov_idx.reshape(-1, 1), bias)
            drift_slope = lr.coef_[0]
        else:
            drift_r, drift_p, drift_slope = np.nan, np.nan, np.nan

        print(f"    MELM ML BIC  = {model_ml.bic:.1f}")
        print(f"    MELM REML var = {model_reml.cov_re.iloc[0, 0]:.6f}")
        print(f"    Drift: r = {drift_r:.4f}, p = {drift_p:.4e}, slope = {drift_slope:.6f}")

        # Collect FOV center coordinates for the spatial schematic
        center_x_col = 'fov_center_x_um' if fov_col == 'fov' else 'fov_perp_center_x_um'
        center_y_col = 'fov_center_y_um' if fov_col == 'fov' else 'fov_perp_center_y_um'

        fov_centers = pd.DataFrame()
        if center_x_col in adata.obs.columns and center_y_col in adata.obs.columns:
            fov_centers = (
                adata.obs[[fov_col, center_x_col, center_y_col]]
                .drop_duplicates()
                .rename(columns={fov_col: 'fov',
                                 center_x_col: 'fov_center_x_um',
                                 center_y_col: 'fov_center_y_um'})
                .sort_values('fov')
            )

        results[label] = {
            'random_intercepts': re_df,
            'drift_r': drift_r,
            'drift_p': drift_p,
            'drift_slope': drift_slope,
            'fov_centers': fov_centers,
            'melm_bic': model_ml.bic,
        }

    return results


# ============================================================================
# Pipeline: run batch effect evaluation across all normalization layers
# Source: 01_run_norm_analyzer.py main loop (lines 267-339)
# ============================================================================

def run_batch_effect_pipeline(
    zarr_path, nextflow_output, technology, tissue_annotation_path,
    output_dir, dataset_name=None, use_log=True, metadata_csv=None,
    n_jobs=None, ci_scheme="cluster",
    he_alignment_path=None, pixel_size=0.2125, n_boot_ci=200,
    vc_column=None, checkpoint_dir=None, fixed_effect_column=None,
):
    """
    Run the full batch effect evaluation across all normalization layers.

    This function orchestrates the complete batch effect analysis:
    1. Load the SpatialData object and FOV identifiers
    2. Load tissue annotations
    3. Compute FOV and tissue summaries from raw data
    4. For each normalization layer:
       a. Load the normalized expression matrix
       b. Fit OLS and MELM (ML for comparison, REML for intercepts)
       c. Compute bootstrap CI for tau-squared
       d. Extract statistics and random intercepts
       e. Unload to free memory
    5. Save all outputs

    Parameters
    ----------
    zarr_path : str
        Path to the QC'd .zarr file.
    nextflow_output : str
        Path to the Nextflow output directory (containing norm/ subdirectory).
    technology : str
        iST platform name ('CosMx', 'Xenium', 'MERSCOPE').
    tissue_annotation_path : str
        Path to tissue annotation file (.csv or .geojson).
    output_dir : str
        Directory to write output files.
    dataset_name : str, optional
        Name for the dataset. If None, inferred from the zarr filename.
    use_log : bool, optional
        Whether to log-transform library sizes before modelling (default True).
    metadata_csv : str, optional
        Path to the enriched metadata CSV (produced by assign_fov.py).
        If provided, FOV identifiers are loaded from this file.
        If None, falls back to searching nextflow_output/norm/ for metadata.
    he_alignment_path : str, optional
        Path to 3x3 H&E-to-Xenium alignment affine CSV.
        Only used for GeoJSON annotations.
    pixel_size : float, optional
        Image pixel size in um (default 0.2125).
    n_boot_ci : int, optional
        Number of bootstrap iterations for tau-squared CI (default 200).
    vc_column : str, optional
        Column name for an additional variance component in the MELM.
    fixed_effect_column : str, optional
        obs column used as the fixed effect instead of the tissue
        annotation (e.g. 'location' for the TMA's primary, metastasis and
        normal cores).
    """
    os.makedirs(output_dir, exist_ok=True)

    if dataset_name is None:
        dataset_name = os.path.basename(zarr_path).replace('.zarr', '')

    # ---- Load data ----
    print(f"\n{'='*60}")
    print(f"BISTRO Batch Effect Evaluation: {dataset_name}")
    print(f"{'='*60}")

    sd_obj = load_spatialdata(zarr_path)

    # Load FOV identifiers from the enriched metadata
    if metadata_csv is not None:
        print(f"Loading FOV from enriched metadata: {metadata_csv}")
        meta = pd.read_csv(metadata_csv, index_col=0)
        adata = sd_obj.tables['filtered']
        meta.index = meta.index.astype(str)
        adata.obs.index = adata.obs.index.astype(str)
        # Load all FOV-related columns (including perpendicular if present)
        fov_cols = [c for c in meta.columns
                    if c.startswith('fov')]
        for col in fov_cols:
            if col in meta.columns:
                adata.obs[col] = meta.loc[adata.obs.index, col].values
        print(f"  FOV loaded: {adata.obs['fov'].nunique()} unique FOVs")
        if 'fov_perp' in adata.obs.columns:
            print(f"  Perpendicular FOV loaded: {adata.obs['fov_perp'].nunique()} unique FOVs")
    else:
        load_fov_from_metadata(sd_obj, nextflow_output)

    os.makedirs(output_dir, exist_ok=True)
    load_tissue_annotations(sd_obj, tissue_annotation_path,
                            he_alignment_path=he_alignment_path,
                            pixel_size_um=pixel_size,
                            polygon_out_path=os.path.join(
                                output_dir,
                                f"{dataset_name}_tissue_polygons.geojson"))

    tissue_column = fixed_effect_column or 'tissue_annotations'
    if tissue_column not in sd_obj.tables['filtered'].obs.columns:
        raise ValueError(f"Fixed-effect column '{tissue_column}' is not in obs")
    print(f"Fixed effect: {tissue_column}")

    # ---- Discover normalization layers ----
    norm_map = discover_norm_layers(nextflow_output)
    if not norm_map:
        print("No normalization layers found. Exiting.")
        return

    # ---- Compute FOV and tissue summaries from raw data ----
    print(f"\n{'='*60}")
    print("Computing FOV and tissue summaries from raw data...")
    print(f"{'='*60}")

    fov_summary, tissue_ls = compute_fov_and_tissue_summaries(
        sd_obj, dataset_name
    )

    fov_path = os.path.join(output_dir, f"{dataset_name}_fov_summary.csv")
    fov_summary.to_csv(fov_path, index=False)
    print(f"Saved FOV summary: {fov_path}")

    if not tissue_ls.empty:
        tissue_path = os.path.join(output_dir, f"{dataset_name}_tissue_ls.csv")
        tissue_ls.to_csv(tissue_path, index=False)
        print(f"Saved tissue LS: {tissue_path}")

    # ---- Save cell overview for the spatial scatter plot in the report ----
    # Subsample to max 50k cells for plotting efficiency
    adata = sd_obj.tables['filtered']
    overview_cols = []
    for xc, yc in [('x_global_um', 'y_global_um'), ('x_local_um', 'y_local_um')]:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            overview_cols = [xc, yc]
            break
    if overview_cols and 'tissue_annotations' in adata.obs.columns:
        # Carry fov and library size alongside the coordinates: the report's
        # diagnostics use them to show FOV assignment and acquisition drift,
        # and to zoom on a single FOV to check that annotation labels form
        # contiguous patches rather than noise.
        extra = [c for c in ['fov'] if c in adata.obs.columns]
        cell_overview = adata.obs[overview_cols + ['tissue_annotations']
                                  + extra].copy()
        cell_overview.columns = ['x', 'y', 'tissue'] + extra

        ls_col = next((c for c in ['nCount_RNA', 'total_counts']
                       if c in adata.obs.columns), None)
        if ls_col:
            cell_overview['library_size'] = adata.obs[ls_col].values
        else:
            try:
                X = adata.X
                cell_overview['library_size'] = (
                    np.asarray(X.sum(axis=1)).ravel())
            except Exception:
                pass
        if len(cell_overview) > 50000:
            cell_overview = cell_overview.sample(n=50000, random_state=42)
        overview_path = os.path.join(output_dir, f"{dataset_name}_cell_overview.csv")
        cell_overview.to_csv(overview_path, index=False)
        print(f"Saved cell overview: {overview_path} ({len(cell_overview)} cells)")

    # ---- Compute drift comparison (row vs perpendicular) on raw data ----
    drift_comparison = compute_drift_comparison(
        sd_obj, dataset_name, technology, use_log=use_log
    )

    if drift_comparison is not None:
        # Save drift comparison results
        drift_rows = []
        all_drift_intercepts = []
        for direction, stats in drift_comparison.items():
            drift_rows.append({
                'dataset': dataset_name,
                'technology': technology,
                'scan_direction': direction,
                'drift_pearson_r': stats['drift_r'],
                'drift_pearson_p': stats['drift_p'],
                'drift_slope': stats['drift_slope'],
                'melm_bic': stats['melm_bic'],
            })
            all_drift_intercepts.append(stats['random_intercepts'])

            if not stats['fov_centers'].empty:
                centers_path = os.path.join(
                    output_dir,
                    f"{dataset_name}_fov_centers_{direction}.csv"
                )
                stats['fov_centers'].to_csv(centers_path, index=False)

        pd.DataFrame(drift_rows).to_csv(
            os.path.join(output_dir, f"{dataset_name}_drift_comparison.csv"),
            index=False
        )
        pd.concat(all_drift_intercepts, ignore_index=True).to_csv(
            os.path.join(output_dir, f"{dataset_name}_drift_intercepts.csv"),
            index=False
        )
        print(f"Saved drift comparison outputs")

    # ---- Set up per-layer checkpoint directory (persistent across restarts) ----
    if checkpoint_dir is None:
        checkpoint_dir = os.path.join(output_dir, 'checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)
    print(f"Per-layer checkpoint directory: {checkpoint_dir}")

    # ---- Process each normalization layer ----
    all_results = {}
    all_random_intercepts = []
    summary_rows = []

    for layer_name in sorted(norm_map.keys()):
        print(f"\n{'='*60}")
        print(f"Processing: {layer_name}")
        print(f"{'='*60}")

        ckpt_summary = os.path.join(checkpoint_dir, f"{dataset_name}_{layer_name}_summary.csv")
        ckpt_re      = os.path.join(checkpoint_dir, f"{dataset_name}_{layer_name}_random_intercepts.csv")
        ckpt_model   = os.path.join(checkpoint_dir, f"{dataset_name}_{layer_name}_model.pkl")

        # Resume: skip layers that already have a complete checkpoint
        if (os.path.exists(ckpt_summary)
                and os.path.exists(ckpt_re)
                and os.path.exists(ckpt_model)):
            print(f"  Checkpoint found for {layer_name}, loading and skipping recomputation.")
            summary_rows.append(pd.read_csv(ckpt_summary).iloc[0].to_dict())
            all_random_intercepts.append(pd.read_csv(ckpt_re))
            with open(ckpt_model, 'rb') as f:
                all_results[layer_name] = pickle.load(f)
            continue

        try:
            # Load layer into memory
            load_layer_into_spatialdata(sd_obj, layer_name, norm_map)

            # Run evaluation
            result = evaluate_batch_effect_single_layer(
                sd_obj, layer_name, dataset_name, technology,
                use_log=use_log, n_boot_ci=n_boot_ci,
                vc_column=vc_column, n_jobs=n_jobs, ci_scheme=ci_scheme,
                tissue_column=tissue_column,
            )

            re_layer = result['random_effects_mixedlm_model_random_intercept']

            summary_row = {
                'layer': result['layer'],
                'dataset': result['dataset'],
                'technology': result['technology'],
                'aic_no_fov': result['aic_no_fov'],
                'bic_no_fov': result['bic_no_fov'],
                'aic_with_fov': result['aic_with_fov'],
                'bic_with_fov': result['bic_with_fov'],
                'aic_melm': result['aic_mixedlm_random_intercept'],
                'bic_melm': result['bic_mixedlm_random_intercept'],
                'delta_bic': result['bic_no_fov'] - result['bic_mixedlm_random_intercept'],
                'lrt_pvalue': result['pvalue_ols_no_fov_vs_mixedlm_random_intercept'],
                'var_random_intercept_ml': result['var_mixedlm_ml'],
                'var_random_intercept_reml': result['var_mixedlm_reml'],
                'var_ci_lower': result['var_ci_lower'],
                'var_ci_upper': result['var_ci_upper'],
                'ci_method': result['ci_method'],
                'drift_pearson_r': result['drift_pearson_r'],
                'drift_pearson_p': result['drift_pearson_p'],
                'drift_slope': result['drift_slope'],
                'drift_slope_se': result.get('drift_slope_se', np.nan),
                'drift_axis': result.get('drift_axis', 'rank'),
                'drift_response': result.get('drift_response', 'relative_bias'),
                'drift_pearson_r_fovid': result.get('drift_pearson_r_fovid', np.nan),
                'drift_pearson_p_fovid': result.get('drift_pearson_p_fovid', np.nan),
                'drift_slope_fovid': result.get('drift_slope_fovid', np.nan),
                'drift_slope_se_fovid': result.get('drift_slope_se_fovid', np.nan),
                'n_nonfinite_intercepts': result.get('n_nonfinite_intercepts', 0),
                'response_is_library_size': result.get('response_is_library_size', True),
                'response_scale': result.get('response_scale', 'log1p'),
                'response_sd': result.get('response_sd', np.nan),
                'response_degenerate': result.get('response_degenerate', False),
                'lrt_reference': result.get('lrt_reference', ''),
                'vc_column': result.get('vc_column', ''),
                'var_vc_reml': result.get('var_vc_reml', np.nan),
                'model_structure': result.get('model_structure', 'fov'),
                'fixed_effect_column': result.get('fixed_effect_column', ''),
            }

            # Write the per-layer checkpoint immediately, before moving on.
            # This is the unit of resumable work: a layer that gets here is
            # never recomputed on a subsequent restart.
            pd.DataFrame([summary_row]).to_csv(ckpt_summary, index=False)
            re_layer.to_csv(ckpt_re, index=False)
            with open(ckpt_model, 'wb') as f:
                pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"  Checkpoint written for {layer_name}")

            all_results[layer_name] = result
            all_random_intercepts.append(re_layer)
            summary_rows.append(summary_row)

            # Unload to free memory
            unload_layer(sd_obj, layer_name)

        except Exception as e:
            print(f"FAILED processing layer {layer_name}: {e}")
            unload_layer(sd_obj, layer_name)
            continue

    # ---- Save outputs ----
    print(f"\n{'='*60}")
    print("Saving outputs...")
    print(f"{'='*60}")

    # 1. Summary CSV
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, f"{dataset_name}_batch_effect_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary: {summary_path}")

    # 2. Random intercepts CSV
    if all_random_intercepts:
        re_df = pd.concat(all_random_intercepts, ignore_index=True)
        re_path = os.path.join(output_dir, f"{dataset_name}_random_intercepts.csv")
        re_df.to_csv(re_path, index=False)
        print(f"Saved random intercepts: {re_path}")

    # 3. Full model objects (pickle) for downstream inspection
    models_path = os.path.join(output_dir, f"{dataset_name}_batch_effect_models.pkl")
    with open(models_path, 'wb') as f:
        pickle.dump(all_results, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved models: {models_path}")

    print(f"\nBatch effect evaluation complete for {dataset_name}.")


# ============================================================================
# CLI
# ============================================================================

def build_parser():
    """Build argument parser for the batch effect evaluation CLI."""
    p = argparse.ArgumentParser(
        description="BISTRO: FOV Batch Effect Evaluation. "
                    "Fits OLS and MELM models per normalization method to "
                    "quantify the FOV-based intra-slide batch effect.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--zarr", required=True,
                   help="Path to the QC'd .zarr file")
    p.add_argument("--nextflow_output", required=True,
                   help="Path to the Nextflow output directory "
                        "(must contain a norm/ subdirectory)")
    p.add_argument("--technology", required=True,
                   choices=["CosMx", "Xenium", "MERSCOPE"],
                   help="iST platform")
    p.add_argument("--tissue_annotation", required=True,
                   help="Path to tissue annotation file (.csv or .geojson)")
    p.add_argument("--output_dir", required=True,
                   help="Directory to write output files")
    p.add_argument("--dataset_name", default=None,
                   help="Dataset name (default: inferred from zarr path)")
    p.add_argument("--metadata", default=None,
                   help="Path to enriched metadata CSV from assign_fov.py "
                        "(if not provided, searches nextflow_output/norm/)")
    p.add_argument("--he_alignment", default=None,
                   help="Path to 3x3 H&E-to-Xenium alignment affine CSV "
                        "(no header). Only used for GeoJSON annotations.")
    p.add_argument("--pixel_size", type=float, default=0.2125,
                   help="Image pixel size in um (default 0.2125)")
    p.add_argument("--n_boot_ci", type=int, default=200,
                   help="Number of bootstrap iterations for tau-squared CI "
                        "(default 200)")
    p.add_argument("--n_jobs", type=int, default=None,
                   help="Parallel workers for the bootstrap. Pass task.cpus; "
                        "the default oversubscribes under SLURM.")
    p.add_argument("--ci_scheme", choices=["cluster", "within"], default="cluster",
                   help="Bootstrap resampling unit for the tau^2 interval. "
                        "'cluster' resamples FOVs (variance-component interval, "
                        "the default); 'within' resamples cells inside a fixed "
                        "FOV set (conditional interval, much narrower).")
    p.add_argument("--use_log", action="store_true", default=True,
                   help="Log-transform library sizes before modelling")
    p.add_argument("--vc_column", default=None,
                   help="Column in obs that FOVs are nested in (e.g. "
                        "'patient' for a TMA). It becomes the grouping factor "
                        "and FOV a variance component within it. If not "
                        "provided, FOV is the only random intercept.")
    p.add_argument("--fixed_effect_column", default=None,
                   help="obs column to use as the fixed effect instead of the "
                        "tissue annotation (e.g. 'location' for a TMA).")
    p.add_argument("--checkpoint_dir", default=None,
                   help="Directory for per-layer checkpoints, persistent "
                        "across job restarts. Layers with a complete "
                        "checkpoint are skipped on resume. Defaults to "
                        "<output_dir>/checkpoints.")
    return p


def main():
    """CLI entry point."""
    args = build_parser().parse_args()
    run_batch_effect_pipeline(
        zarr_path=args.zarr,
        nextflow_output=args.nextflow_output,
        technology=args.technology,
        tissue_annotation_path=args.tissue_annotation,
        output_dir=args.output_dir,
        dataset_name=args.dataset_name,
        use_log=args.use_log,
        metadata_csv=args.metadata,
        he_alignment_path=args.he_alignment,
        pixel_size=args.pixel_size,
        n_boot_ci=args.n_boot_ci,
        vc_column=args.vc_column,
        checkpoint_dir=args.checkpoint_dir,
        fixed_effect_column=args.fixed_effect_column,
        n_jobs=args.n_jobs,
        ci_scheme=args.ci_scheme,
    )


if __name__ == "__main__":
    main()