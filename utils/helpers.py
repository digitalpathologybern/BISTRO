"""
BISTRO utilities: shared helper functions and constants.

All functions in this module are pure (stateless) and shared across
the evaluation scripts. Sources are noted per function.

Coordinate convention
---------------------
All spatial coordinates are in **micrometers (um)**. Column names:
  - x_local_um, y_local_um   : cell centroid within its FOV
  - x_global_um, y_global_um : cell centroid in slide coordinates
  - fov_center_x_um, fov_center_y_um : FOV center in slide coordinates
  - area_um2                  : cell area in square micrometers
"""

import numpy as np
import pandas as pd
from scipy.stats import chi2


# ============================================================================
# Constants
# Source: 02_transformation.ipynb, cell 48
# ============================================================================

TRANSFORM_ORDER = [
    "raw", "log(x+1)", "log(x+0.01)", "log(x+0.1)",
    "log(x+0.5)", "log(x+10)", "sqrt",
    "acosh", "pearson_residuals",
]

NORM_ORDER = [
    "none", "cpm", "cp100", "cp10k", "scran",
    "deseq2", "tmm", "areaNorm",
    "spanorm.logpac", "spanorm.pearson", "scTransform",
]

RAW_ONLY_NORMS = {"scTransform", "spanorm.logpac", "spanorm.pearson"}
COMPOUND_NORMS = {"spanorm-logpac", "spanorm-pearson"}


# ============================================================================
# Library size computation
# Source: NormAnalyzerPipeline._calculate_ls (norm_analyzer.py, line 244)
# ============================================================================

def calculate_library_size(matrix):
    """
    Compute per-cell library size (sum of counts across genes).

    Handles both dense numpy arrays and scipy sparse matrices.

    Parameters
    ----------
    matrix : np.ndarray or scipy.sparse matrix
        Cell-by-gene expression matrix.

    Returns
    -------
    np.ndarray
        1-D array of library sizes, one per cell.
    """
    if hasattr(matrix, 'A1') or 'scipy' in str(type(matrix)):
        library_sizes = matrix.sum(axis=1).A1
    else:
        library_sizes = matrix.sum(axis=1)

    return library_sizes


# ============================================================================
# Statistical helpers for batch effect evaluation
# Source: nested functions inside NormAnalyzerPipeline._process_layer_fast
#         (norm_analyzer.py, lines 1694-1705)
# ============================================================================

def var_ci_chisq(var_hat, n_groups, alpha=0.05):
    """
    Chi-squared approximation for confidence interval of random intercept
    variance (tau-squared).

    Implements Equation 5 from the BISTRO manuscript.

    Parameters
    ----------
    var_hat : float
        Estimated variance of the random intercept.
    n_groups : int
        Number of FOV groups.
    alpha : float
        Significance level (default 0.05 for 95% CI).

    Returns
    -------
    np.ndarray
        Array [lower_bound, upper_bound] of the confidence interval.
    """
    df = max(n_groups - 1, 1)
    lower = (df * var_hat) / chi2.ppf(1 - alpha / 2, df)
    upper = (df * var_hat) / chi2.ppf(alpha / 2, df)
    return np.array([lower, upper])


def _bootstrap_refit_worker(args):
    """
    One bootstrap iteration for the random-intercept variance.

    scheme='cluster' resamples FOVs with replacement, relabelling each draw with
    a fresh synthetic group id. scheme='within' resamples cells inside a fixed
    FOV set.

    Parameters
    ----------
    args : tuple
        (boot_idx, df, formula, group_col, seed, scheme, vc_formula)

    Returns
    -------
    float
        tau-squared from the refitted model, or np.nan on failure.
    """
    import statsmodels.formula.api as smf

    boot_idx, df, formula, group_col, seed, scheme, vc_formula = args
    rng = np.random.default_rng(seed + boot_idx)

    if scheme == "cluster":
        groups = list(df.groupby(group_col).indices.items())
        picks = rng.choice(len(groups), size=len(groups), replace=True)
        parts = []
        for new_id, gi in enumerate(picks):
            part = df.iloc[groups[gi][1]].copy()
            part["_boot_group"] = new_id
            parts.append(part)
        df_boot = pd.concat(parts, ignore_index=True)
        fit_groups = df_boot["_boot_group"]
    else:
        parts = []
        for _, group_df in df.groupby(group_col):
            n = len(group_df)
            parts.append(group_df.iloc[rng.choice(n, size=n, replace=True)])
        df_boot = pd.concat(parts, ignore_index=True)
        fit_groups = df_boot[group_col]

    try:
        kwargs = dict(formula=formula, groups=fit_groups, data=df_boot,
                      re_formula="~1")
        if vc_formula:
            kwargs["vc_formula"] = vc_formula
        model = smf.mixedlm(**kwargs).fit(method=["powell"], reml=True,
                                          maxiter=200)
        return model.cov_re.iloc[0, 0]
    except Exception:
        return np.nan


def bootstrap_var_ci(df, formula, group_col, n_boot=200, alpha=0.05,
                     seed=42, n_jobs=None, scheme="cluster", vc_formula=None):
    """
    Parametric bootstrap confidence interval for the random intercept
    variance (tau-squared) from a mixed-effects linear model.

    scheme selects the resampling unit. vc_formula, when given, is carried into
    every refit.

    This replaces the chi-squared approximation (var_ci_chisq) which has
    poor coverage for small group counts.

    Parameters
    ----------
    df : pd.DataFrame
        The modelling dataframe with dependent variable, fixed effects,
        and group column.
    formula : str
        The model formula (e.g. 'log_LS ~ C(tissue_annotations)').
    group_col : str
        Column name for the grouping variable (FOV).
    n_boot : int
        Number of bootstrap iterations (default 200).
    alpha : float
        Significance level (default 0.05 for 95% CI).
    seed : int
        Random seed for reproducibility.
    n_jobs : int or None
        Number of parallel workers. Pass task.cpus under SLURM.
    scheme : {'cluster', 'within'}
        Resampling unit. 'cluster' resamples FOVs, 'within' resamples cells.
    vc_formula : dict or None
        Extra variance components, passed through to every refit.

    Returns
    -------
    np.ndarray
        Array [lower_bound, upper_bound] of the bootstrap CI.
        Returns [np.nan, np.nan] if bootstrap fails.
    """
    import os
    from concurrent.futures import ProcessPoolExecutor

    if n_jobs is None:
        n_jobs = min(n_boot, os.cpu_count() or 4)

    # Prepare worker arguments
    worker_args = [
        (b, df, formula, group_col, seed, scheme, vc_formula)
        for b in range(n_boot)
    ]

    # Run bootstrap in parallel
    tau_sq_samples = []
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_bootstrap_refit_worker, worker_args)
        tau_sq_samples = [r for r in results if not np.isnan(r)]

    if len(tau_sq_samples) < 10:
        print(f"  WARNING: only {len(tau_sq_samples)} successful bootstrap "
              f"iterations (out of {n_boot}), CI may be unreliable")
        if len(tau_sq_samples) == 0:
            return np.array([np.nan, np.nan])

    tau_sq_arr = np.array(tau_sq_samples)
    lower = np.percentile(tau_sq_arr, 100 * alpha / 2)
    upper = np.percentile(tau_sq_arr, 100 * (1 - alpha / 2))

    return np.array([lower, upper])


def lrt_pvalue(model_restricted, model_full, n_vc=1, boundary=True):
    """
    Likelihood ratio test p-value comparing a restricted model to a full model.

    Implements Equation 3 from the BISTRO manuscript:
        Lambda = 2 * (loglik_full - loglik_restricted)

    Referenced to the boundary-corrected mixture
        0.5 * chi2_{n_vc} + 0.5 * chi2_{n_vc - 1}
    (Self and Liang 1987; Stram and Lee 1994).

    Parameters
    ----------
    model_restricted : statsmodels results object
        The restricted (null) model (e.g., OLS without FOV).
    model_full : statsmodels results object
        The full model (e.g., MELM with FOV random intercept).
    n_vc : int
        Number of variance components the full model adds over the restricted one.
    boundary : bool
        Use the boundary-corrected mixture reference. Set False to reproduce
        the naive chi2_1 test.

    Returns
    -------
    float
        p-value from the likelihood ratio test.
    """
    lr_stat = 2 * (model_full.llf - model_restricted.llf)
    lr_stat = max(lr_stat, 0.0)
    if not boundary:
        return chi2.sf(lr_stat, n_vc)
    upper = chi2.sf(lr_stat, n_vc)
    lower = chi2.sf(lr_stat, n_vc - 1) if n_vc > 1 else 0.0
    return 0.5 * (upper + lower)


# ============================================================================
# FOV size definitions (in micrometers)
# Source: NormAnalyzerPipeline._get_fov_size (norm_analyzer.py, line 211)
#
# Updated: values are now in micrometers (um), matching the standardized
# coordinate schema where all spatial columns are in um.
# ============================================================================

# Xenium rasterisation pitch in um. Override per dataset with params.fovTileUm.
XENIUM_TILE_UM = [600.0, 720.0]
XENIUM_TILE_UM_PROTOTYPE = [600.0, 875.0]


def get_fov_size(technology, scale=1, tile_um=None):
    """
    Return the FOV tile dimensions [width_um, height_um] for a given
    iST technology, used for pseudo-FOV rasterization.

    All values are in **micrometers**.


    Parameters
    ----------
    technology : str
        One of 'CosMx', 'Xenium', 'MERSCOPE'.
    scale : float
        Scaling factor applied to the FOV dimensions (default 1).
    tile_um : sequence of two floats, optional
        Explicit [width_um, height_um] pitch, overriding the default.

    Returns
    -------
    list
        [width_um, height_um].
    """
    if tile_um is not None:
        w, h = float(tile_um[0]), float(tile_um[1])
        return [w / scale, h / scale]
    if technology == 'Xenium':
        return [XENIUM_TILE_UM[0] / scale, XENIUM_TILE_UM[1] / scale]
    elif technology == 'CosMx':
        return [510.72, 510.72]
    elif technology == 'MERSCOPE':
        # TODO: verify MERSCOPE pixel size and derive FOV dimensions in um
        return [223 / scale, 223 / scale]
    else:
        raise ValueError(f"Unknown technology: {technology}")


# ============================================================================
# FOV rasterization for Xenium / MERSCOPE
# Source: rasterize_gene_expression standalone function (norm_analyzer.py, line 16)
#
# Updated: coordinates are now in um (no /1000 conversion), tile sizes in um,
# output center columns use new naming convention.
# ============================================================================

def rasterize_fov_grid(meta, tile_size_um, min_cells_per_fov=0,
                       scan_axis='row', fov_col='fov'):
    """
    Assign pseudo-FOV identifiers to cells by binning spatial coordinates
    into a regular grid. Used for platforms that do not natively report
    FOV indices (Xenium, MERSCOPE).

    The scan_axis parameter controls the direction of FOV numbering:
      - 'row' (default): iterate y-bins top-to-bottom, within each row
        iterate x-bins left-to-right. Horizontal raster scan.
      - 'col' (perpendicular): iterate x-bins left-to-right, within each
        column iterate y-bins top-to-bottom. Vertical raster scan.

    If min_cells_per_fov > 0, small FOVs are merged with the preceding
    FOV (used only for SpaNorm compatibility). Default (0): every tile
    gets its own FOV.

    Parameters
    ----------
    meta : pd.DataFrame
        Cell metadata with columns 'x_local_um' and 'y_local_um'
        (coordinates in micrometers).
    tile_size_um : list
        [width_um, height_um] of each pseudo-FOV tile, in micrometers.
    min_cells_per_fov : int
        Minimum cells per FOV; smaller FOVs are merged (default 0 = no merging).
    scan_axis : str
        Rasterization direction: 'row' (y-primary, default) or 'col'
        (x-primary, perpendicular).
    fov_col : str
        Name of the output FOV column (default 'fov'). Use e.g. 'fov_perp'
        to store a perpendicular assignment alongside the default.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added columns:
          '{fov_col}', 'fov_center_x_um', 'fov_center_y_um'
        (or '{fov_col}_center_x_um', '{fov_col}_center_y_um' if fov_col != 'fov').
    """
    meta = meta.copy()

    # Bin coordinates directly in um (no unit conversion needed)
    meta['_x_bin'] = np.floor(meta['x_local_um'] / tile_size_um[0]).astype(int)
    meta['_y_bin'] = np.floor(meta['y_local_um'] / tile_size_um[1]).astype(int)

    # Build list of FOV keys in raster order.
    fov_keys_ordered = []

    if scan_axis == 'row':
        # Row-primary: iterate y (top-to-bottom), then x (left-to-right)
        unique_y_bins = sorted(meta['_y_bin'].unique())
        for y in unique_y_bins:
            x_bins = sorted(meta.loc[meta['_y_bin'] == y, '_x_bin'].unique())
            for x in x_bins:
                fov_keys_ordered.append(f"{x}_{y}")
    elif scan_axis == 'col':
        # Column-primary (perpendicular): iterate x (left-to-right), then y (top-to-bottom)
        unique_x_bins = sorted(meta['_x_bin'].unique())
        for x in unique_x_bins:
            y_bins = sorted(meta.loc[meta['_x_bin'] == x, '_y_bin'].unique())
            for y in y_bins:
                fov_keys_ordered.append(f"{x}_{y}")
    else:
        raise ValueError(f"scan_axis must be 'row' or 'col', got '{scan_axis}'")

    # Assign temporary fov ids
    temp_fov_ids = {key: idx for idx, key in enumerate(fov_keys_ordered)}
    meta['_key'] = meta['_x_bin'].astype(str) + '_' + meta['_y_bin'].astype(str)
    meta['_temp_fov'] = meta['_key'].map(temp_fov_ids)

    # Count number of cells per temporary FOV
    fov_counts = (
        meta.groupby('_temp_fov')
        .size()
        .reset_index(name='cell_count')
        .sort_values('_temp_fov')
    )

    # Assign FOV labels. If min_cells_per_fov > 0, small FOVs are merged
    # into the preceding FOV (used only for SpaNorm compatibility).
    # Default (min_cells_per_fov=0): every tile gets its own FOV.
    fov_label_map = {}
    fov_counter = 0
    prev_valid_fov = 0

    for i, row in fov_counts.iterrows():
        this_fov = row['_temp_fov']
        count = row['cell_count']
        if min_cells_per_fov == 0 or fov_counter == 0 or count >= min_cells_per_fov:
            fov_label_map[this_fov] = fov_counter
            prev_valid_fov = fov_counter
            fov_counter += 1
        else:
            fov_label_map[this_fov] = prev_valid_fov

    # Assign final FOV labels and centers using the caller-specified column names
    center_x_col = 'fov_center_x_um' if fov_col == 'fov' else f'{fov_col}_center_x_um'
    center_y_col = 'fov_center_y_um' if fov_col == 'fov' else f'{fov_col}_center_y_um'

    meta[fov_col] = meta['_temp_fov'].map(fov_label_map)
    meta[center_x_col] = (meta['_x_bin'] + 0.5) * tile_size_um[0]
    meta[center_y_col] = (meta['_y_bin'] + 0.5) * tile_size_um[1]

    # Cleanup temporary columns
    meta.drop(columns=['_x_bin', '_y_bin', '_key', '_temp_fov'], inplace=True)

    return meta


# ============================================================================
# Layer name parsing
# Source: 02_transformation.ipynb, cell 48 (line 171)
# ============================================================================

def parse_layer(layer):
    """
    Split a compound layer name into (normalization, transformation).

    Handles special compound normalization names like 'spanorm-logpac'
    which should be parsed as normalization='spanorm.logpac', transformation='raw',
    rather than splitting at the first hyphen.

    Parameters
    ----------
    layer : str
        Layer name, e.g. 'cpm-log(x+1)', 'spanorm-logpac', 'scTransform'.

    Returns
    -------
    tuple
        (normalization, transformation) strings.
    """
    layer = str(layer)
    for cn in sorted(COMPOUND_NORMS, key=len, reverse=True):
        if layer == cn:
            return cn.replace("-", "."), "raw"
        if layer.startswith(cn + "-"):
            return cn.replace("-", "."), layer[len(cn) + 1:]
    if "-" not in layer:
        return layer, "raw"
    norm, trans = layer.split("-", 1)
    return norm, trans