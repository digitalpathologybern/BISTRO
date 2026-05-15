#!/usr/bin/env python
"""
BISTRO Report Generator
========================

Collects all evaluation outputs and produces a self-contained HTML report.

Usage
-----
    python generate_report.py \\
        --results_dir /path/to/evaluation_output/ \\
        --dataset_name MyDataset \\
        --output_html /path/to/BISTRO_report.html
"""

import argparse
import base64
import io
import os
import sys
from glob import glob
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'utils'))
from helpers import NORM_ORDER, TRANSFORM_ORDER


# ============================================================================
# Utility
# ============================================================================

def fig_to_base64(fig, dpi=150):
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return encoded


def find_file(results_dir, pattern):
    """Find a single file matching a glob pattern (searches subdirs too)."""
    matches = glob(os.path.join(results_dir, pattern))
    if matches:
        return matches[0]
    matches = glob(os.path.join(results_dir, '**', pattern), recursive=True)
    return matches[0] if matches else None


# ============================================================================
# Section 0: Cell Overview
# ============================================================================

def plot_cell_overview(overview_df):
    """Scatter plot of cells colored by tissue annotation."""
    if overview_df is None or overview_df.empty:
        return None
    if 'x' not in overview_df.columns or 'y' not in overview_df.columns:
        return None
    if 'tissue' not in overview_df.columns:
        return None

    overview_df = overview_df.dropna(subset=['tissue'])
    overview_df['tissue'] = overview_df['tissue'].astype(str)

    # Remove unannotated cells
    overview_df = overview_df[overview_df['tissue'] != 'None']
    if overview_df.empty:
        return None

    tissues = sorted(overview_df['tissue'].unique())
    cmap = plt.colormaps['tab20']
    color_map = {t: cmap(i % 20) for i, t in enumerate(tissues)}

    fig, ax = plt.subplots(figsize=(8, 7))

    for tissue in tissues:
        sub = overview_df[overview_df['tissue'] == tissue]
        ax.scatter(sub['x'], sub['y'], s=0.3, alpha=0.5,
                   color=color_map[tissue], label=tissue, rasterized=True)

    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Cell Overview (colored by tissue annotation)', fontsize=11)

    # Legend outside plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5),
              fontsize=7, markerscale=10, frameon=False)

    plt.tight_layout()
    return fig


# ============================================================================
# Section 1: Batch Effect figures
# ============================================================================

def plot_fov_heatmap(fov_df):
    """Rasterized 2D heatmap of mean library size per FOV."""
    if fov_df is None or fov_df.empty:
        return None
    if 'fov_center_x_um' not in fov_df.columns or 'fov_center_y_um' not in fov_df.columns:
        return None

    fig, ax = plt.subplots(figsize=(7, 6))
    norm = mcolors.Normalize(vmin=fov_df['mean_ls'].min(), vmax=fov_df['mean_ls'].max())
    cmap = plt.cm.plasma

    # Estimate FOV tile size from spacing
    xs = sorted(fov_df['fov_center_x_um'].unique())
    ys = sorted(fov_df['fov_center_y_um'].unique())
    dx = np.median(np.diff(xs)) if len(xs) > 1 else 0.5
    dy = np.median(np.diff(ys)) if len(ys) > 1 else 0.5

    for _, row in fov_df.iterrows():
        x = row['fov_center_x_um'] - dx / 2
        y = row['fov_center_y_um'] - dy / 2
        color = cmap(norm(row['mean_ls']))
        rect = Rectangle((x, y), dx, dy, linewidth=0.3,
                          edgecolor='white', facecolor=color)
        ax.add_patch(rect)

    ax.set_xlim(fov_df['fov_center_x_um'].min() - dx,
                fov_df['fov_center_x_um'].max() + dx)
    ax.set_ylim(fov_df['fov_center_y_um'].min() - dy,
                fov_df['fov_center_y_um'].max() + dy)
    ax.set_aspect('equal')
    ax.axis('off')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02, label='Mean Library Size')
    ax.set_title('Mean Library Size per FOV')
    plt.tight_layout()
    return fig


def plot_tissue_barplot(tissue_df):
    """Vertical barplot of mean library size per tissue region."""
    if tissue_df is None or tissue_df.empty:
        return None

    df = tissue_df.sort_values('mean_ls', ascending=False)
    df['tissue_region'] = df['tissue_region'].astype(str)
    fig, ax = plt.subplots(figsize=(max(4, len(df) * 0.8), 5))
    bars = ax.bar(df['tissue_region'], df['mean_ls'],
                  color='#4C72B0', edgecolor='black', linewidth=0.5)

    if 'std_ls' in df.columns:
        ax.errorbar(df['tissue_region'], df['mean_ls'], yerr=df['std_ls'],
                     fmt='none', ecolor='black', capsize=3, linewidth=0.8)

    ax.set_ylabel('Mean Library Size')
    ax.set_title('Mean Library Size per Tissue Region')
    ax.tick_params(axis='x', rotation=45)
    for label in ax.get_xticklabels():
        label.set_ha('right')
    plt.tight_layout()
    return fig


def plot_delta_bic(summary_df):
    """Bar plot of delta-BIC per normalization with significance markers."""
    if summary_df is None or summary_df.empty:
        return None

    df = summary_df.sort_values('delta_bic', ascending=True)
    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.4)))
    colors = ['#2ca02c' if v > 0 else '#d62728' for v in df['delta_bic']]
    ax.barh(df['layer'], df['delta_bic'], color=colors, edgecolor='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Delta-BIC (positive = MELM preferred)')
    ax.set_title('FOV Batch Effect: Model Comparison (OLS vs MELM)')

    for i, (_, row) in enumerate(df.iterrows()):
        qval = row.get('lrt_qvalue', row.get('lrt_pvalue', 1.0))
        if qval < 0.001:
            marker = '***'
        elif qval < 0.01:
            marker = '**'
        elif qval < 0.05:
            marker = '*'
        else:
            marker = 'n.s.'
        x_pos = row['delta_bic'] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.01
        ax.text(x_pos, i, marker, va='center', fontsize=8, color='red')

    plt.tight_layout()
    return fig


def plot_drift_and_sd(intercepts_df):
    """
    Side-by-side: (left) random intercept drift scatter with r, p, slope;
    (right) barplot of SD of random intercepts for normalization 'none'.
    """
    if intercepts_df is None or intercepts_df.empty:
        return None

    layers = intercepts_df['layer'].unique()
    target = 'none' if 'none' in layers else layers[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5),
                                    gridspec_kw={'width_ratios': [2.5, 1]})

    # ---- Left: drift scatter for the target layer ----
    df_target = intercepts_df[intercepts_df['layer'] == target].sort_values('fov')
    fov_idx = np.arange(len(df_target))
    bias = df_target['relative_bias'].values

    ax1.scatter(fov_idx, bias, s=20, alpha=0.7, edgecolors='black', linewidth=0.3)

    if len(fov_idx) > 2 and np.std(bias) > 0:
        from scipy.stats import pearsonr
        r_val, p_val = pearsonr(fov_idx, bias)
        z = np.polyfit(fov_idx, bias, 1)
        p_line = np.poly1d(z)
        ax1.plot(fov_idx, p_line(fov_idx), 'r-', linewidth=1.5)
        ax1.text(0.02, 0.95,
                 f'slope = {z[0]:.4f}\nr = {r_val:.3f}\np = {p_val:.2e}',
                 transform=ax1.transAxes, fontsize=9, va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax1.axhline(0, color='grey', linestyle='--', linewidth=0.5)
    ax1.set_xlabel('FOV Acquisition Index')
    ax1.set_ylabel('Relative Bias (exp(u_j) - 1)')
    ax1.set_title(f'FOV Signal Drift (layer: {target})')

    # ---- Right: SD barplot across normalizations for "none" ----
    sd_per_layer = (
        intercepts_df.groupby('layer')['random_intercept']
        .std()
        .sort_values(ascending=False)
        .reset_index()
        .rename(columns={'random_intercept': 'sd'})
    )

    colors = ['#d62728' if l == target else '#4C72B0' for l in sd_per_layer['layer']]
    ax2.barh(sd_per_layer['layer'], sd_per_layer['sd'],
             color=colors, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('SD of Random Intercepts')
    ax2.set_title('Random Intercept Variability')
    ax2.invert_yaxis()

    plt.tight_layout()
    return fig


def _draw_raster_schematic(ax, fov_centers, scan_direction, title):
    """
    Draw a schematic of the FOV grid with arrows showing the rasterization
    direction and FOV numbering order.
    """
    if fov_centers is None or fov_centers.empty:
        ax.text(0.5, 0.5, 'No FOV coordinates', ha='center', va='center',
                transform=ax.transAxes, fontsize=9, color='grey')
        ax.set_title(title)
        return

    centers = fov_centers.sort_values('fov')
    n_fovs = len(centers)

    # Color by FOV index
    norm = mcolors.Normalize(vmin=0, vmax=max(n_fovs - 1, 1))
    cmap = plt.cm.viridis

    # Estimate tile size from spacing
    xs = sorted(centers['fov_center_x_um'].unique())
    ys = sorted(centers['fov_center_y_um'].unique())
    dx = np.median(np.diff(xs)) if len(xs) > 1 else 0.5
    dy = np.median(np.diff(ys)) if len(ys) > 1 else 0.5

    for _, row in centers.iterrows():
        x = row['fov_center_x_um'] - dx / 2
        y = row['fov_center_y_um'] - dy / 2
        color = cmap(norm(row['fov']))
        rect = Rectangle((x, y), dx, dy, linewidth=0.3,
                          edgecolor='white', facecolor=color)
        ax.add_patch(rect)

    # Draw arrows showing scan direction.
    # Prioritize row/column transitions (jumps between rows or columns),
    # which are the most informative for understanding the scan pattern.
    # Always show the first transition and at least one more.
    if n_fovs > 1:
        # Classify each consecutive pair as "inline" or "transition"
        transitions = []  # indices where a row/col jump occurs
        inlines = []      # indices within the same row/col
        for i in range(n_fovs - 1):
            r0 = centers.iloc[i]
            r1 = centers.iloc[i + 1]
            ddx = abs(r1['fov_center_x_um'] - r0['fov_center_x_um'])
            ddy = abs(r1['fov_center_y_um'] - r0['fov_center_y_um'])
            # A transition is when movement occurs in both axes or
            # primarily in the non-scan axis
            if scan_direction == 'row':
                is_transition = ddy > dy * 0.5  # significant Y jump = new row
            else:
                is_transition = ddx > dx * 0.5  # significant X jump = new column
            if is_transition:
                transitions.append(i)
            else:
                inlines.append(i)

        # Always draw the first transition and subsample the rest
        arrows_to_draw = set()
        if transitions:
            arrows_to_draw.add(transitions[0])  # first row/col transition
            if len(transitions) > 1:
                arrows_to_draw.add(transitions[1])  # second transition
            # Add a few more transitions spread across the slide
            step = max(1, len(transitions) // 4)
            for idx in transitions[::step]:
                arrows_to_draw.add(idx)

        # Also draw a few inline arrows to show within-row direction
        inline_step = max(1, len(inlines) // 6)
        for idx in inlines[::inline_step]:
            arrows_to_draw.add(idx)

        for i in sorted(arrows_to_draw):
            r0 = centers.iloc[i]
            r1 = centers.iloc[i + 1]
            ax.annotate('', xy=(r1['fov_center_x_um'], r1['fov_center_y_um']),
                        xytext=(r0['fov_center_x_um'], r0['fov_center_y_um']),
                        arrowprops=dict(arrowstyle='->', color='red',
                                        lw=1.0, alpha=0.7))

    ax.set_xlim(centers['fov_center_x_um'].min() - dx,
                centers['fov_center_x_um'].max() + dx)
    ax.set_ylim(centers['fov_center_y_um'].min() - dy,
                centers['fov_center_y_um'].max() + dy)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=10, fontweight='bold')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02, label='FOV index')


def _draw_drift_scatter(ax, intercepts_df, drift_r, drift_p, drift_slope, title):
    """Draw a drift scatter plot on the given axes."""
    df = intercepts_df.sort_values('fov')
    fov_idx = np.arange(len(df))
    bias = df['relative_bias'].values

    ax.scatter(fov_idx, bias, s=15, alpha=0.7, edgecolors='black', linewidth=0.3)

    if len(fov_idx) > 2 and np.isfinite(drift_slope):
        z = np.polyfit(fov_idx, bias, 1)
        p_line = np.poly1d(z)
        ax.plot(fov_idx, p_line(fov_idx), 'r-', linewidth=1.5)
        ax.text(0.02, 0.95,
                f'slope = {drift_slope:.4f}\nr = {drift_r:.3f}\np = {drift_p:.2e}',
                transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)
    ax.set_xlabel('FOV Acquisition Index', fontsize=8)
    ax.set_ylabel('Relative Bias', fontsize=8)
    ax.set_title(title, fontsize=10, fontweight='bold')


def _draw_sd_barplot(ax, intercepts_df, title):
    """Draw horizontal barplot of SD of random intercepts per normalization."""
    sd_per_layer = (
        intercepts_df.groupby('layer')['random_intercept']
        .std()
        .sort_values(ascending=False)
        .reset_index()
        .rename(columns={'random_intercept': 'sd'})
    )
    if sd_per_layer.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title(title)
        return

    colors = ['#4C72B0'] * len(sd_per_layer)
    # Highlight 'none' normalization
    for i, l in enumerate(sd_per_layer['layer']):
        if l == 'none':
            colors[i] = '#d62728'

    ax.barh(sd_per_layer['layer'], sd_per_layer['sd'],
            color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xlabel('SD of Random Intercepts', fontsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.invert_yaxis()
    ax.tick_params(axis='y', labelsize=6)


def plot_drift_comparison(drift_comp_df, drift_intercepts_df,
                          intercepts_per_norm_df,
                          fov_centers_row=None, fov_centers_col=None):
    """
    Two-row figure comparing row-primary and column-primary FOV orderings.

    Each row contains 3 columns:
      - Left: rasterization schematic (FOV grid colored by acquisition order
        with arrows showing scan direction)
      - Center: drift scatter plot (relative bias vs FOV index with r, p, slope)
      - Right: SD of random intercepts per normalization for this FOV ordering

    Only generated for rasterized platforms (Xenium, MERSCOPE).

    Parameters
    ----------
    drift_comp_df : pd.DataFrame
        Drift comparison summary with columns: scan_direction, drift_pearson_r,
        drift_pearson_p, drift_slope.
    drift_intercepts_df : pd.DataFrame
        Per-FOV random intercepts for both orderings, with column
        'scan_direction' ('row' or 'col').
    intercepts_per_norm_df : pd.DataFrame
        Per-FOV random intercepts from the main batch effect analysis
        (all normalization layers), used for the SD barplot. Must have
        'layer' and 'random_intercept' columns.
    fov_centers_row : pd.DataFrame or None
        FOV center coordinates for row-primary ordering.
    fov_centers_col : pd.DataFrame or None
        FOV center coordinates for col-primary ordering.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    if drift_comp_df is None or drift_comp_df.empty:
        return None
    if drift_intercepts_df is None or drift_intercepts_df.empty:
        return None

    fig, axes = plt.subplots(2, 3, figsize=(18, 10),
                              gridspec_kw={'width_ratios': [1, 1.8, 0.8]})

    for i, (direction, label, centers) in enumerate([
        ('row', 'Row-primary scan (horizontal)', fov_centers_row),
        ('col', 'Column-primary scan (vertical)', fov_centers_col),
    ]):
        row_data = drift_comp_df[drift_comp_df['scan_direction'] == direction]
        re_data = drift_intercepts_df[
            drift_intercepts_df['scan_direction'] == direction
        ]

        if row_data.empty or re_data.empty:
            for j in range(3):
                axes[i, j].text(0.5, 0.5, f'No data for {direction}',
                                ha='center', va='center',
                                transform=axes[i, j].transAxes)
            continue

        r_val = row_data.iloc[0]['drift_pearson_r']
        p_val = row_data.iloc[0]['drift_pearson_p']
        slope = row_data.iloc[0]['drift_slope']

        # Left: rasterization schematic
        _draw_raster_schematic(axes[i, 0], centers, direction,
                               f'{label}\n(FOV ordering)')

        # Center: drift scatter
        _draw_drift_scatter(axes[i, 1], re_data, r_val, p_val, slope,
                            f'{label}\n(Drift: r={r_val:.3f}, p={p_val:.2e})')

        # Right: SD barplot for this direction
        if intercepts_per_norm_df is not None and not intercepts_per_norm_df.empty:
            _draw_sd_barplot(axes[i, 2], intercepts_per_norm_df,
                             f'Intercept Variability\n({direction}-primary)')
        else:
            axes[i, 2].text(0.5, 0.5, 'No normalization data',
                            ha='center', va='center',
                            transform=axes[i, 2].transAxes)

    fig.suptitle('FOV Signal Drift & Random Intercept Variability',
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()
    return fig


# ============================================================================
# Section 2: Transformation figures
# ============================================================================

def plot_transformation_heatmap(summary_df, metric='binned_cv'):
    """Heatmap of a variance stabilization metric across norm x transform."""
    if summary_df is None or summary_df.empty:
        return None

    pivot = summary_df.pivot_table(
        index='normalization', columns='transformation',
        values=metric, aggfunc='median')

    row_order = [r for r in NORM_ORDER if r in pivot.index]
    col_order = [c for c in TRANSFORM_ORDER if c in pivot.columns]
    pivot = pivot.reindex(index=row_order, columns=col_order)

    fig, ax = plt.subplots(figsize=(max(6, len(col_order) * 0.9),
                                     max(4, len(row_order) * 0.5)))
    im = ax.imshow(pivot.values, aspect='auto', cmap='RdYlGn_r')
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order, fontsize=7)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=6)

    plt.colorbar(im, ax=ax, shrink=0.8, label=metric)
    ax.set_title(f'Variance Stabilization: {metric}')
    plt.tight_layout()
    return fig


def plot_pc1_heatmap(corr_df):
    """Heatmap of PC1-mean correlation across norm x transform."""
    if corr_df is None or corr_df.empty:
        return None

    pivot = corr_df.pivot_table(
        index='normalization', columns='transformation',
        values='correlation', aggfunc='median')

    row_order = [r for r in NORM_ORDER if r in pivot.index]
    col_order = [c for c in TRANSFORM_ORDER if c in pivot.columns]
    pivot = pivot.reindex(index=row_order, columns=col_order)

    fig, ax = plt.subplots(figsize=(max(6, len(col_order) * 0.9),
                                     max(4, len(row_order) * 0.5)))
    im = ax.imshow(pivot.values, aspect='auto', cmap='RdYlGn_r',
                   vmin=0, vmax=1)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order, fontsize=7)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=6)

    plt.colorbar(im, ax=ax, shrink=0.8, label='Spearman r(mean, |PC1|)')
    ax.set_title('PC1-Mean Expression Correlation')
    plt.tight_layout()
    return fig


# ============================================================================
# Section 3: HVG Benchmark figures
# ============================================================================

def plot_hvg_vs_random(hvg_df):
    """
    HVG vs random gene set comparison (Figure 4C style).

    For each normalization method (rows) and HVG fraction (columns),
    shows a violin of the random gene set ARI values with the HVG ARI
    overlaid as a colored diamond.

    Parameters
    ----------
    hvg_df : pd.DataFrame
        HVG benchmark results with columns: layer, hvg_fraction,
        hvg_stability, random_stability_values (JSON string of individual
        random ARI scores).

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    import json as _json

    if hvg_df is None or hvg_df.empty:
        return None
    if 'hvg_stability' not in hvg_df.columns or 'random_stability_values' not in hvg_df.columns:
        return None

    # Filter to rows that have both HVG and random data, exclude frac=1.0
    df = hvg_df.dropna(subset=['hvg_stability', 'random_stability_values']).copy()
    df = df[df['hvg_fraction'] < 1.0]
    if df.empty:
        return None

    # Parse JSON random values
    def _parse_random(val):
        if isinstance(val, str):
            try:
                return _json.loads(val)
            except Exception:
                return []
        return []

    df['random_values'] = df['random_stability_values'].apply(_parse_random)

    layers = sorted(df['layer'].unique())
    fractions = sorted(df['hvg_fraction'].unique())

    n_layers = len(layers)
    n_fracs = len(fractions)

    if n_layers == 0 or n_fracs == 0:
        return None

    # Create subplot grid: one row per normalization, one column per HVG fraction
    fig, axes = plt.subplots(
        n_layers, n_fracs,
        figsize=(max(3 * n_fracs, 6), max(2 * n_layers, 4)),
        squeeze=False, sharey=True
    )

    cmap = plt.colormaps['tab10']

    for i, layer in enumerate(layers):
        for j, frac in enumerate(fractions):
            ax = axes[i, j]
            row = df[(df['layer'] == layer) & (df['hvg_fraction'] == frac)]

            if row.empty:
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                        transform=ax.transAxes, fontsize=8, color='grey')
            else:
                row = row.iloc[0]
                rand_vals = row['random_values']
                hvg_val = row['hvg_stability']

                # Violin of random values
                if len(rand_vals) > 1:
                    parts = ax.violinplot(rand_vals, positions=[0], showmeans=False,
                                          showmedians=True, widths=0.6)
                    for pc in parts['bodies']:
                        pc.set_facecolor('#BBBBBB')
                        pc.set_alpha(0.7)
                    for key in ['cmins', 'cmaxes', 'cmedians', 'cbars']:
                        if key in parts:
                            parts[key].set_color('grey')
                elif len(rand_vals) == 1:
                    ax.scatter([0], rand_vals, color='#BBBBBB', s=30, zorder=2)

                # HVG diamond
                ax.scatter([0], [hvg_val], marker='D', color=cmap(i),
                           s=60, zorder=3, edgecolors='black', linewidth=0.5)

            ax.set_xticks([])

            # Column title (top row only)
            if i == 0:
                ax.set_title(f'{frac:.0%} HVG', fontsize=8)

            # Row label (first column only)
            if j == 0:
                ax.set_ylabel(layer, fontsize=7, rotation=0, ha='right',
                              va='center', labelpad=50)

    fig.suptitle('HVG vs Random Gene Set: Clustering Stability (ARI)',
                 fontsize=11, fontweight='bold', y=1.02)

    # Add a legend explaining the symbols
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#4C72B0',
               markeredgecolor='black', markersize=8, label='HVG'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#BBBBBB',
               markersize=8, label='Random gene sets'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2,
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout()
    return fig


def plot_hvg_benchmark_panels(hvg_df):
    """
    Side-by-side: (left) spatial coherence vs HVG fraction,
    (right) mean AUROC vs HVG fraction, grouped by normalization.
    """
    if hvg_df is None or hvg_df.empty:
        return None

    # Filter to needed columns
    needed_cols = {'hvg_fraction', 'layer'}
    has_coherence = 'spatial_coherence' in hvg_df.columns
    has_auroc = 'marker_auroc_mean' in hvg_df.columns

    if not has_coherence and not has_auroc:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Get unique layers for coloring
    layers = sorted(hvg_df['layer'].unique())
    cmap = plt.colormaps['tab10']
    color_map = {l: cmap(i) for i, l in enumerate(layers)}

    # ---- Left: Spatial Coherence ----
    ax = axes[0]
    if has_coherence:
        for layer in layers:
            sub = hvg_df[hvg_df['layer'] == layer].sort_values('hvg_fraction')
            sub = sub.dropna(subset=['spatial_coherence'])
            if not sub.empty:
                ax.plot(sub['hvg_fraction'], sub['spatial_coherence'],
                        'o-', color=color_map[layer], label=layer,
                        markersize=5, linewidth=1.2)
        ax.set_xlabel('HVG Fraction')
        ax.set_ylabel('Spatial Coherence')
        ax.set_title('Clustering Spatial Coherence')
        ax.legend(fontsize=6, ncol=2, loc='best')
    else:
        ax.text(0.5, 0.5, 'No spatial coherence data', ha='center',
                va='center', transform=ax.transAxes)
        ax.set_title('Clustering Spatial Coherence')

    # ---- Right: Marker AUROC ----
    ax = axes[1]
    if has_auroc:
        for layer in layers:
            sub = hvg_df[hvg_df['layer'] == layer].sort_values('hvg_fraction')
            sub = sub.dropna(subset=['marker_auroc_mean'])
            if not sub.empty:
                ax.plot(sub['hvg_fraction'], sub['marker_auroc_mean'],
                        's-', color=color_map[layer], label=layer,
                        markersize=5, linewidth=1.2)
        ax.set_xlabel('HVG Fraction')
        ax.set_ylabel('Mean AUROC')
        ax.set_title('Mean AUROC of Marker Gene Discriminability')
        ax.legend(fontsize=6, ncol=2, loc='best')
    else:
        ax.text(0.5, 0.5, 'No AUROC data', ha='center',
                va='center', transform=ax.transAxes)
        ax.set_title('Mean AUROC of Marker Gene Discriminability')

    plt.tight_layout()
    return fig


# ============================================================================
# Section 4: Annotation agreement figure
# ============================================================================

def plot_ari_heatmap(ari_df):
    """Heatmap of pairwise ARI between normalization methods."""
    if ari_df is None or ari_df.empty:
        return None

    fig, ax = plt.subplots(figsize=(max(6, len(ari_df.columns) * 0.7),
                                     max(5, len(ari_df.index) * 0.6)))
    im = ax.imshow(ari_df.values.astype(float), aspect='auto',
                   cmap='YlOrRd', vmin=0, vmax=1)
    ax.set_xticks(range(len(ari_df.columns)))
    ax.set_xticklabels(ari_df.columns, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(ari_df.index)))
    ax.set_yticklabels(ari_df.index, fontsize=7)

    for i in range(ari_df.shape[0]):
        for j in range(ari_df.shape[1]):
            val = float(ari_df.values[i, j])
            if np.isfinite(val):
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=6)

    plt.colorbar(im, ax=ax, shrink=0.8, label='ARI')
    ax.set_title('Cell Type Annotation Agreement (Pairwise ARI)')
    plt.tight_layout()
    return fig


# ============================================================================
# HTML generation
# ============================================================================

def generate_html_report(results_dir, dataset_name, output_html):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sections = []

    # ==== Section 0: Cell Overview ====
    cell_overview_file = find_file(results_dir, "*cell_overview*")
    if cell_overview_file:
        cell_overview = pd.read_csv(cell_overview_file)
        fig = plot_cell_overview(cell_overview)
        if fig:
            b64 = fig_to_base64(fig)
            sections.append(
                f'<h2>Spatial Overview</h2>\n'
                f'<img src="data:image/png;base64,{b64}" />\n'
            )

    # ==== Section 1: Batch Effect ====
    batch_summary_file = find_file(results_dir, "*batch_effect_summary*")
    intercepts_file = find_file(results_dir, "*random_intercepts*")
    fov_summary_file = find_file(results_dir, "*fov_summary*")
    tissue_ls_file = find_file(results_dir, "*tissue_ls*")

    batch_summary = pd.read_csv(batch_summary_file) if batch_summary_file else None
    
    # Apply Benjamini-Hochberg FDR correction across all LRT p-values
    if batch_summary is not None and 'lrt_pvalue' in batch_summary.columns:
        from statsmodels.stats.multitest import multipletests
        valid_mask = batch_summary['lrt_pvalue'].notna()
        if valid_mask.sum() > 1:
            _, qvals, _, _ = multipletests(
                batch_summary.loc[valid_mask, 'lrt_pvalue'].values,
                method='fdr_bh'
            )
            batch_summary.loc[valid_mask, 'lrt_qvalue'] = qvals
        else:
            batch_summary['lrt_qvalue'] = batch_summary['lrt_pvalue']
    
    intercepts = pd.read_csv(intercepts_file) if intercepts_file else None
    fov_summary = pd.read_csv(fov_summary_file) if fov_summary_file else None
    tissue_ls = pd.read_csv(tissue_ls_file) if tissue_ls_file else None

    section_html = '<h2>1. FOV Batch Effect Evaluation</h2>\n'

    # 1a. FOV heatmap and tissue barplot side by side
    fig_fov = plot_fov_heatmap(fov_summary)
    fig_tissue = plot_tissue_barplot(tissue_ls)
    if fig_fov or fig_tissue:
        section_html += '<div style="display:flex; gap:20px; flex-wrap:wrap; align-items:flex-start;">\n'
        if fig_fov:
            b64 = fig_to_base64(fig_fov)
            section_html += f'<div><h3>Mean Library Size per FOV</h3><img src="data:image/png;base64,{b64}" /></div>\n'
        if fig_tissue:
            b64 = fig_to_base64(fig_tissue)
            section_html += f'<div><h3>Mean Library Size per Tissue Region</h3><img src="data:image/png;base64,{b64}" /></div>\n'
        section_html += '</div>\n'

    # 1b. Delta-BIC
    if batch_summary is not None:
        fig = plot_delta_bic(batch_summary)
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<h3>Delta-BIC: OLS vs MELM</h3>\n'
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'

        section_html += '<h3>Summary Table</h3>\n'
        section_html += batch_summary.round(4).to_html(
            index=False, classes='dataframe', border=0)

    # 1c. Drift analysis
    drift_comp_file = find_file(results_dir, "*drift_comparison*")
    drift_intercepts_file = find_file(results_dir, "*drift_intercepts*")
    fov_centers_row_file = find_file(results_dir, "*fov_centers_row*")
    fov_centers_col_file = find_file(results_dir, "*fov_centers_col*")

    if drift_comp_file and drift_intercepts_file:
        # Rasterized platform: merged 2x3 figure (schematic | drift | SD per direction)
        drift_comp = pd.read_csv(drift_comp_file)
        drift_int = pd.read_csv(drift_intercepts_file)
        centers_row = pd.read_csv(fov_centers_row_file) if fov_centers_row_file else None
        centers_col = pd.read_csv(fov_centers_col_file) if fov_centers_col_file else None

        fig = plot_drift_comparison(drift_comp, drift_int, intercepts,
                                    centers_row, centers_col)
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<h3>FOV Signal Drift & Random Intercept Variability</h3>\n'
            section_html += ('<p style="color:#555; font-size:12px;">'
                             'Because the true microscope acquisition order is unknown '
                             'for this platform, both row-primary and column-primary '
                             'rasterization orderings are evaluated. The ordering that '
                             'reveals a stronger drift (higher |r|) is more likely to '
                             'align with the actual scan direction.</p>\n')
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'
    else:
        # Native FOV platform (CosMx): show standard drift + SD plot
        if intercepts is not None:
            fig = plot_drift_and_sd(intercepts)
            if fig:
                b64 = fig_to_base64(fig)
                section_html += f'<h3>FOV Signal Drift & Random Intercept Variability</h3>\n'
                section_html += f'<img src="data:image/png;base64,{b64}" />\n'

    sections.append(section_html)

    # ==== Section 2: Transformation Analysis ====
    trans_summary_file = find_file(results_dir, "*transformation_summary*")
    pc1_corr_file = find_file(results_dir, "*pc1_correlation*")

    trans_summary = pd.read_csv(trans_summary_file) if trans_summary_file else None
    pc1_corr = pd.read_csv(pc1_corr_file) if pc1_corr_file else None

    section_html = '<h2>2. Variance Stabilization Analysis</h2>\n'

    # Only binned_cv heatmap (removed slope, spearman, lowess_r2)
    if trans_summary is not None:
        fig = plot_transformation_heatmap(trans_summary, metric='binned_cv')
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<h3>Binned CV of Variance</h3>\n'
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'

    if pc1_corr is not None:
        fig = plot_pc1_heatmap(pc1_corr)
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<h3>PC1-Mean Correlation</h3>\n'
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'

    sections.append(section_html)

    # ==== Section 3: HVG Benchmark ====
    hvg_files = glob(os.path.join(results_dir, '**', '*.csv'), recursive=True)
    hvg_files = [f for f in hvg_files if 'checkpoint' in f.lower() or 'benchmark' in f.lower()]

    section_html = '<h2>3. HVG Selection Benchmark</h2>\n'

    hvg_combined = None
    if hvg_files:
        hvg_dfs = []
        for f in hvg_files:
            try:
                hvg_dfs.append(pd.read_csv(f))
            except Exception:
                continue
        if hvg_dfs:
            hvg_combined = pd.concat(hvg_dfs, ignore_index=True)

    if hvg_combined is not None and not hvg_combined.empty:
        # HVG vs random gene set comparison (Figure 4C style)
        fig = plot_hvg_vs_random(hvg_combined)
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<h3>HVG vs Random Gene Set: Clustering Stability</h3>\n'
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'

        # Spatial coherence + AUROC side-by-side
        fig = plot_hvg_benchmark_panels(hvg_combined)
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'

        # Summary table
        cols_to_show = [c for c in [
            'layer', 'hvg_fraction', 'clustering_stability_median_ari',
            'hvg_stability', 'random_stability_mean',
            'spatial_coherence', 'marker_auroc_mean',
        ] if c in hvg_combined.columns]

        if cols_to_show:
            section_html += '<h3>Results Summary</h3>\n'
            section_html += hvg_combined[cols_to_show].round(4).to_html(
                index=False, classes='dataframe', border=0)
    else:
        section_html += '<p>No HVG benchmark results found.</p>\n'

    sections.append(section_html)

    # ==== Section 4: Annotation Agreement ====
    ari_file = find_file(results_dir, "*pairwise_ari*")

    section_html = '<h2>4. Cell Type Annotation Agreement</h2>\n'
    if ari_file:
        ari_df = pd.read_csv(ari_file, index_col=0)
        fig = plot_ari_heatmap(ari_df)
        if fig:
            b64 = fig_to_base64(fig)
            section_html += f'<img src="data:image/png;base64,{b64}" />\n'
        section_html += '<h3>Pairwise ARI Matrix</h3>\n'
        section_html += ari_df.round(4).to_html(classes='dataframe', border=0)
    else:
        section_html += '<p>No annotation agreement results found.</p>\n'

    sections.append(section_html)

    # ==== Assemble HTML ====
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>BISTRO Report: {dataset_name}</title>
    <style>
        body {{
            font-family: Arial, Helvetica, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #fafafa;
            color: #333;
            line-height: 1.5;
        }}
        h1 {{
            border-bottom: 3px solid #2c3e50;
            padding-bottom: 10px;
            color: #2c3e50;
        }}
        h2 {{
            color: #2c3e50;
            border-bottom: 1px solid #bdc3c7;
            padding-bottom: 5px;
            margin-top: 40px;
        }}
        h3 {{ color: #555; }}
        img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            margin: 10px 0;
        }}
        table.dataframe {{
            border-collapse: collapse;
            width: 100%;
            font-size: 12px;
            margin: 15px 0;
        }}
        table.dataframe th, table.dataframe td {{
            border: 1px solid #ddd;
            padding: 6px 10px;
            text-align: right;
        }}
        table.dataframe th {{
            background-color: #2c3e50;
            color: white;
        }}
        table.dataframe tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        .meta {{
            color: #888;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <h1>BISTRO Report: {dataset_name}</h1>
    <p class="meta">Generated: {timestamp}</p>
    <p class="meta">
        BISTRO (Bias Identification in Spatial Transcriptomics) evaluates
        FOV batch effects, variance stabilization, HVG selection, and
        cell type annotation consistency across normalization methods.
    </p>
    <hr>
    {''.join(sections)}
    <hr>
    <p class="meta">End of BISTRO report.</p>
</body>
</html>"""

    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)
    with open(output_html, 'w') as f:
        f.write(html)
    print(f"Report saved: {output_html}")


# ============================================================================
# CLI
# ============================================================================

def build_parser():
    p = argparse.ArgumentParser(description="BISTRO: Generate HTML Report",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--results_dir", required=True)
    p.add_argument("--dataset_name", default="BISTRO")
    p.add_argument("--output_html", required=True)
    return p

def main():
    args = build_parser().parse_args()
    generate_html_report(
        results_dir=args.results_dir,
        dataset_name=args.dataset_name,
        output_html=args.output_html)

if __name__ == "__main__":
    main()