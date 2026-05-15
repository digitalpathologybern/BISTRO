#!/usr/bin/env python
"""
BISTRO Annotation Agreement
============================

Computes pairwise Adjusted Rand Index (ARI) between cell type annotations
produced by InSituType under different normalization methods.

Implements the analysis described in Section 2.4 (Figure 4D) of the
BISTRO manuscript.

Outputs
-------
  {dataset}_pairwise_ari_{hvg_threshold}.csv
      Square matrix of ARI values (normalization x normalization).

Usage
-----
    python annotation_agreement.py \\
        --annotation_dir /path/to/nextflow_output/annotation/ \\
        --hvg_threshold 0.75 \\
        --output_dir /path/to/output/
"""

import argparse
import os
import sys
from glob import glob
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score


def load_annotations(annotation_dir, hvg_threshold):
    """
    Load all InSituType annotation CSVs for a given HVG threshold.

    Expects filenames matching the pattern:
        top_{threshold}_HVG_{method}_annotation.csv

    Parameters
    ----------
    annotation_dir : str
        Directory containing annotation CSV files.
    hvg_threshold : str
        HVG threshold string as it appears in filenames (e.g. '0.75').

    Returns
    -------
    dict
        Mapping of {normalization_method: pd.Series of cell type labels}.
    """
    pattern = os.path.join(annotation_dir, f"*{hvg_threshold}*_annotation.csv")
    files = sorted(glob(pattern))

    if not files:
        print(f"Warning: no annotation files found matching: {pattern}")
        return {}

    annotations = {}
    for fpath in files:
        basename = os.path.basename(fpath)
        # Extract method name from filename pattern:
        # e.g. "top_0.75_HVG_cpm_annotation.csv" -> "cpm"
        parts = basename.replace("_annotation.csv", "").split("_HVG_")
        if len(parts) == 2:
            method = parts[1]
        else:
            # Fallback: try to extract from filename
            method = basename.replace("_annotation.csv", "").split("_")[-1]

        df = pd.read_csv(fpath, index_col=0)

        # InSituType output has column 'sup.clust'
        if 'sup.clust' in df.columns:
            annotations[method] = df['sup.clust']
        elif len(df.columns) == 1:
            annotations[method] = df.iloc[:, 0]
        else:
            print(f"  Warning: unexpected columns in {basename}: {df.columns.tolist()}")
            continue

        print(f"  Loaded: {method} ({len(annotations[method])} cells)")

    return annotations


def compute_pairwise_ari(annotations):
    """
    Compute pairwise ARI between all normalization methods.

    For each pair, only cells present in both annotation sets are compared.

    Parameters
    ----------
    annotations : dict
        Mapping of {method_name: pd.Series of cell type labels}.

    Returns
    -------
    pd.DataFrame
        Square DataFrame with ARI values, indexed and columned by method name.
    """
    methods = sorted(annotations.keys())
    n = len(methods)

    ari_matrix = pd.DataFrame(
        np.eye(n), index=methods, columns=methods
    )

    for m1, m2 in combinations(methods, 2):
        s1 = annotations[m1]
        s2 = annotations[m2]

        # Align on shared cell indices
        shared_idx = s1.index.intersection(s2.index)
        if len(shared_idx) < 10:
            print(f"  Warning: only {len(shared_idx)} shared cells for "
                  f"{m1} vs {m2}, setting ARI=NaN")
            ari_matrix.loc[m1, m2] = np.nan
            ari_matrix.loc[m2, m1] = np.nan
            continue

        ari = adjusted_rand_score(
            s1.loc[shared_idx].values,
            s2.loc[shared_idx].values,
        )
        ari_matrix.loc[m1, m2] = round(ari, 4)
        ari_matrix.loc[m2, m1] = round(ari, 4)

    return ari_matrix


def run_annotation_agreement(annotation_dir, hvg_threshold, output_dir,
                             dataset_name=None):
    """
    Load annotations for a given HVG threshold and compute pairwise ARI.

    Parameters
    ----------
    annotation_dir : str
        Directory containing InSituType annotation CSVs.
    hvg_threshold : str
        HVG threshold (e.g. '0.75').
    output_dir : str
        Directory to write output.
    dataset_name : str, optional
        Dataset name for output file naming.
    """
    os.makedirs(output_dir, exist_ok=True)

    if dataset_name is None:
        dataset_name = "dataset"

    print(f"Loading annotations for HVG threshold: {hvg_threshold}")
    annotations = load_annotations(annotation_dir, hvg_threshold)

    if len(annotations) < 2:
        print("Need at least 2 normalization methods. Exiting.")
        return

    print(f"\nComputing pairwise ARI across {len(annotations)} methods...")
    ari_matrix = compute_pairwise_ari(annotations)

    output_path = os.path.join(
        output_dir,
        f"{dataset_name}_pairwise_ari_{hvg_threshold}.csv"
    )
    ari_matrix.to_csv(output_path)
    print(f"\nSaved: {output_path}")
    print(f"\nARI matrix:\n{ari_matrix.to_string()}")


# ============================================================================
# CLI
# ============================================================================

def build_parser():
    p = argparse.ArgumentParser(
        description="BISTRO: Cell Type Annotation Agreement (Pairwise ARI)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--annotation_dir", required=True,
                   help="Directory containing InSituType annotation CSVs")
    p.add_argument("--hvg_threshold", type=str, default="0.75",
                   help="HVG threshold to evaluate (as in filename)")
    p.add_argument("--output_dir", required=True,
                   help="Directory to write output")
    p.add_argument("--dataset_name", default=None,
                   help="Dataset name for output file naming")
    return p


def main():
    args = build_parser().parse_args()
    run_annotation_agreement(
        annotation_dir=args.annotation_dir,
        hvg_threshold=args.hvg_threshold,
        output_dir=args.output_dir,
        dataset_name=args.dataset_name,
    )


if __name__ == "__main__":
    main()