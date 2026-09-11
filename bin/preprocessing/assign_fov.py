#!/usr/bin/env python
"""
BISTRO FOV Assignment
=====================

Centralized FOV identifier assignment that runs once, immediately
after readZarr. All downstream processes (SpaNorm, batch effect
evaluation, etc.) read FOV from the enriched metadata instead of
running their own rasterization.

Logic:
  1. Read the metadata CSV output by readZarr.
  2. Check if the 'fov' column exists with meaningful values (>1 unique).
  3. If yes (e.g. CosMx): keep native FOVs, compute fov_center_x_um/y_um
     if missing.
  4. If no (Xenium, MERSCOPE): rasterize pseudo-FOVs using the technology-
     specific tile size (in um). Generates both row-primary ('fov') and
     col-primary ('fov_perp') orderings.
  5. Write the enriched metadata CSV and a JSON info file.

Outputs
-------
  {output_dir}/{metadata_basename}  (enriched with fov, fov_center_x_um, fov_center_y_um)
  {output_dir}/fov_info.json        (FOV source info for the report)

Coordinate convention
---------------------
All spatial coordinates are in micrometers (um).
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'utils'))
from helpers import (rasterize_fov_grid, get_fov_size,
                     validate_coordinate_scale)


def check_fov_exists(meta):
    """
    Check whether the metadata contains valid FOV identifiers.

    A valid FOV column must:
      - Exist as a column named 'fov'
      - Have at least some non-null values
      - Contain more than 1 unique value (a single value is likely a placeholder)

    Parameters
    ----------
    meta : pd.DataFrame
        Cell metadata.

    Returns
    -------
    bool
        True if valid FOV identifiers are present.
    """
    if 'fov' not in meta.columns:
        return False
    if meta['fov'].isna().all():
        return False
    if meta['fov'].nunique() <= 1:
        return False
    return True


def compute_fov_centers_from_cells(meta, fov_column='fov'):
    """
    Compute FOV center coordinates as the mean cell position within each FOV.

    Used for platforms with native FOV identifiers (e.g. CosMx) that
    may not include precomputed FOV center coordinates.

    Parameters
    ----------
    meta : pd.DataFrame
        Cell metadata with FOV identifiers and spatial coordinates (um).
    fov_column : str
        Column name for the FOV identifier.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with 'fov_center_x_um' and 'fov_center_y_um' added.
    """
    # Determine which coordinate columns are available (um-unit preferred)
    coord_pairs = [
        ('x_global_um', 'y_global_um'),
        ('x_local_um', 'y_local_um'),
        ('x_centroid', 'y_centroid'),
    ]

    x_col, y_col = None, None
    for xc, yc in coord_pairs:
        if xc in meta.columns and yc in meta.columns:
            x_col, y_col = xc, yc
            break

    if x_col is None:
        print("  Warning: no coordinate columns found, cannot compute FOV centers")
        return meta

    # Compute mean position per FOV
    fov_centers = (
        meta.groupby(fov_column, observed=True)
        .agg(fov_center_x_um=(x_col, 'mean'), fov_center_y_um=(y_col, 'mean'))
    )

    # Assign by mapping on the FOV column rather than merging. A pandas merge
    # discards the index, which would replace the cell IDs with a RangeIndex
    # and break every downstream stage that aligns this metadata against
    # adata.obs.index.
    for col in ['fov_center_x_um', 'fov_center_y_um']:
        meta[col] = meta[fov_column].map(fov_centers[col])

    return meta


def assign_fov(metadata_path, technology, output_dir, fov_map_path=None,
               fov_tile_um=None, check_coord_scale=True):
    """
    Main FOV assignment logic.

    Reads the metadata CSV, checks for existing FOV identifiers, and
    either keeps native FOVs or generates pseudo-FOVs via rasterization.
    Writes the enriched metadata CSV and a JSON info file.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata CSV from readZarr.
    technology : str
        iST platform name ('CosMx', 'Xenium', 'MERSCOPE').
    output_dir : str
        Directory to write the enriched metadata and fov_info.json.

    Returns
    -------
    dict
        FOV info dictionary (also written to fov_info.json).
    """
    print(f"Reading metadata: {metadata_path}")
    meta = pd.read_csv(metadata_path, index_col=0)
    print(f"  Cells: {len(meta)}")
    print(f"  Columns: {list(meta.columns)}")

    coord_check = validate_coordinate_scale(
        meta, technology, raise_on_fail=check_coord_scale)
    if coord_check.get('median_nn_um') is not None:
        print(f"  Median nearest-neighbour distance: "
              f"{coord_check['median_nn_um']:.2f} um "
              f"({coord_check['coord_columns'][0]}/"
              f"{coord_check['coord_columns'][1]})")

    if fov_map_path:
        # ---- Vendor-native FOV from a per-cell map, highest precedence ----
        fmap = pd.read_csv(fov_map_path, index_col=0)
        if 'fov' not in fmap.columns:
            raise ValueError(f"{fov_map_path} has no 'fov' column")
        matched = meta.index.intersection(fmap.index)
        cover = len(matched) / max(len(meta), 1)
        print(f"  FOV map: {fov_map_path}")
        print(f"  Matched {len(matched)} of {len(meta)} cells ({100 * cover:.2f}%)")

        join_keys = pd.Index(meta.index.astype(str))
        if cover < 0.5:
            stripped = pd.Index(
                meta.index.astype(str).str.replace(r'_\d+$', '', regex=True))
            alt = stripped.intersection(fmap.index)
            alt_cover = len(alt) / max(len(meta), 1)
            print(f"  Retry without the region suffix: matched {len(alt)} "
                  f"({100 * alt_cover:.2f}%)")
            if alt_cover > cover:
                join_keys = stripped
                matched, cover = alt, alt_cover
        if cover < 0.5:
            raise ValueError(
                f"FOV map matches only {100 * cover:.2f}% of cells. The cell ids "
                f"do not line up; check that the map was built from the same "
                f"vendor bundle as this zarr.")
        meta['fov'] = fmap['fov'].reindex(join_keys).to_numpy()
        if 'fov_name' in fmap.columns:
            meta['fov_name'] = fmap['fov_name'].reindex(join_keys).to_numpy()
        unmatched = int(meta['fov'].isna().sum())
        if unmatched:
            print(f"  WARNING: {unmatched} cells absent from the map; "
                  f"they carry no FOV and will be dropped downstream")
        meta = compute_fov_centers_from_cells(meta)
        n_fovs = int(meta['fov'].nunique())

        has_perp = 'fov_perp' in fmap.columns
        n_fovs_perp = None
        if has_perp:
            meta['fov_perp'] = fmap['fov_perp'].reindex(join_keys).to_numpy()
            meta = compute_fov_centers_from_cells(meta, fov_column='fov_perp')
            meta = meta.rename(columns={
                'fov_center_x_um': 'fov_perp_center_x_um',
                'fov_center_y_um': 'fov_perp_center_y_um'})
            meta = compute_fov_centers_from_cells(meta, fov_column='fov')
            n_fovs_perp = int(meta['fov_perp'].nunique())

        fov_info = {
            'fov_source': 'vendor_native',
            'technology': technology,
            'n_fovs': n_fovs,
            'n_fovs_perp': n_fovs_perp,
            'fov_map': str(fov_map_path),
            'cells_matched': int(len(matched)),
            'cells_unmatched': unmatched,
            'has_perpendicular': bool(has_perp),
            'message': (f'Vendor-native FOV identities taken from {fov_map_path} '
                        f'({n_fovs} tiles). The integer order is the ASSUMED '
                        f'row-major scan, the same assumption the rasteriser '
                        f'makes; no acquisition timestamp exists.')
        }

    elif check_fov_exists(meta):
        # ---- Native FOV (e.g. CosMx) ----
        n_fovs = int(meta['fov'].nunique())
        print(f"  Native FOV identifiers found: {n_fovs} FOVs")

        # Compute FOV centers if not already present
        if 'fov_center_x_um' not in meta.columns or 'fov_center_y_um' not in meta.columns:
            print("  Computing FOV center coordinates from cell positions...")
            meta = compute_fov_centers_from_cells(meta)

        fov_info = {
            'fov_source': ('rasterized_inherited'
                       if 'fov_perp' in meta.columns else 'native'),
            'technology': technology,
            'n_fovs': n_fovs,
            'message': f'Native FOV identifiers detected ({n_fovs} FOVs). '
                       f'Technology: {technology}.'
        }

    else:
        # ---- Rasterization needed (Xenium, MERSCOPE) ----
        fov_size = get_fov_size(technology, tile_um=fov_tile_um)
        print(f"  No valid FOV identifiers found. Rasterizing pseudo-FOVs...")
        print(f"  Technology: {technology}, tile size: {fov_size} um")

        # Default rasterization: row-primary (y first, then x)
        print(f"  Rasterizing with row-primary scan (default)...")
        meta = rasterize_fov_grid(meta, fov_size, scan_axis='row', fov_col='fov')
        n_fovs = int(meta['fov'].nunique())
        print(f"    {n_fovs} FOVs assigned (row-primary)")

        # Perpendicular rasterization: col-primary (x first, then y)
        print(f"  Rasterizing with col-primary scan (perpendicular)...")
        meta = rasterize_fov_grid(meta, fov_size, scan_axis='col', fov_col='fov_perp')
        n_fovs_perp = int(meta['fov_perp'].nunique())
        print(f"    {n_fovs_perp} FOVs assigned (col-primary)")

        fov_info = {
            'fov_source': 'rasterized',
            'technology': technology,
            'fov_tile_size_um': fov_size,
            'n_fovs': n_fovs,
            'n_fovs_perp': n_fovs_perp,
            'has_perpendicular': True,
            'message': f'Pseudo-FOV identifiers generated by rasterization '
                       f'({n_fovs} FOVs row-primary, {n_fovs_perp} FOVs col-primary, '
                       f'tile size {fov_size} um). Technology: {technology}.'
        }

    # ---- Write enriched metadata to output directory ----
    output_metadata_path = os.path.join(output_dir, os.path.basename(metadata_path))
    meta.to_csv(output_metadata_path)
    print(f"  Enriched metadata written: {output_metadata_path}")
    print(f"  FOV columns: fov={meta['fov'].nunique()} unique, "
          f"fov_center_x_um={'yes' if 'fov_center_x_um' in meta.columns else 'no'}, "
          f"fov_center_y_um={'yes' if 'fov_center_y_um' in meta.columns else 'no'}")

    fov_info['median_nn_um'] = coord_check.get('median_nn_um')
    fov_info['coord_scale_ok'] = coord_check.get('ok')

    # ---- Write FOV info JSON for the report ----
    os.makedirs(output_dir, exist_ok=True)
    info_path = os.path.join(output_dir, 'fov_info.json')
    with open(info_path, 'w') as f:
        json.dump(fov_info, f, indent=2)
    print(f"  FOV info written: {info_path}")
    print(f"  {fov_info['message']}")

    return fov_info


# ============================================================================
# CLI
# ============================================================================

def build_parser():
    """Build argument parser for the FOV assignment CLI."""
    p = argparse.ArgumentParser(
        description="BISTRO: Centralized FOV Assignment. "
                    "Checks for existing FOV identifiers in the metadata and "
                    "generates pseudo-FOVs via rasterization if needed.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--metadata", required=True,
                   help="Path to the metadata CSV from readZarr")
    p.add_argument("--technology", required=True,
                   choices=["CosMx", "Xenium", "MERSCOPE"],
                   help="iST platform")
    p.add_argument("--fov_tile_um", default=None,
                   help="Rasterisation pitch as 'width,height' in um, "
                        "overriding the per-technology default.")
    p.add_argument("--fov_map", default=None,
                   help="Per-cell native FOV map CSV (cell_id, fov, fov_name) "
                        "from bin/preprocessing/make_fov_map.py. Takes "
                        "precedence over the obs fov column and rasterisation.")
    p.add_argument("--output_dir", required=True,
                   help="Directory to write enriched metadata and fov_info.json")
    p.add_argument("--no_coord_check", action="store_true",
                   help="Warn instead of failing when the coordinate scale "
                        "check finds a median nearest-neighbour distance "
                        "outside the plausible micrometer range.")
    return p


def main():
    """CLI entry point."""
    args = build_parser().parse_args()
    assign_fov(
        metadata_path=args.metadata,
        technology=args.technology,
        output_dir=args.output_dir,
        fov_map_path=args.fov_map,
        fov_tile_um=([float(v) for v in args.fov_tile_um.split(',')]
                     if args.fov_tile_um else None),
        check_coord_scale=not args.no_coord_check,
    )


if __name__ == "__main__":
    main()