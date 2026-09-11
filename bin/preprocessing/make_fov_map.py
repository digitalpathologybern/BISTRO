#!/usr/bin/env python
"""
Build a per-cell native FOV map from a vendor transcripts table.

Xenium reports its tile identity per transcript as `fov_name` (analysis software
1.4 onward), not in the cell table. This assigns each cell its modal tile and
emits an integer `fov` ordered row-major over measured tile positions.

Outputs:
  <out>.csv    cell_id, fov, fov_name
  <out>.json   fov_name -> fov mapping, ordering rule, straddling fraction

Usage:
  python make_fov_map.py --transcripts transcripts.parquet --output fov_map
"""
import argparse
import json
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

UNASSIGNED = {"UNASSIGNED", "-1", ""}


def _decode(v):
    return v.decode() if isinstance(v, (bytes, bytearray)) else v


def build_map(transcripts_path, verbose=True):
    import pyarrow.parquet as pq

    pf = pq.ParquetFile(transcripts_path)
    cols = pf.schema_arrow.names
    if "fov_name" not in cols:
        raise SystemExit(
            f"{transcripts_path} has no 'fov_name' column. This vendor export "
            f"predates it (Xenium analysis software 1.4). Columns: {cols}"
        )

    pair_counts = Counter()          # (cell_id, fov_name) -> n transcripts
    tile_x = defaultdict(list)       # fov_name -> sampled x
    tile_y = defaultdict(list)

    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=["cell_id", "fov_name",
                                           "x_location", "y_location"])
        df = pd.DataFrame({
            "cell_id": [_decode(v) for v in t.column("cell_id").to_pylist()],
            "fov_name": [_decode(v) for v in t.column("fov_name").to_pylist()],
            "x": t.column("x_location").to_pylist(),
            "y": t.column("y_location").to_pylist(),
        })
        df = df[~df["cell_id"].isin(UNASSIGNED)]
        if df.empty:
            continue
        for (c, f), k in df.groupby(["cell_id", "fov_name"]).size().items():
            pair_counts[(c, f)] += int(k)
        g = df.groupby("fov_name")[["x", "y"]].min()
        for f, row in g.iterrows():
            tile_x[f].append(row["x"])
            tile_y[f].append(row["y"])
        if verbose:
            print(f"  row group {rg + 1}/{pf.num_row_groups}", flush=True)

    # Modal tile per cell.
    best, total, straddle = {}, defaultdict(int), 0
    per_cell = defaultdict(list)
    for (c, f), k in pair_counts.items():
        per_cell[c].append((k, f))
        total[c] += k
    for c, lst in per_cell.items():
        lst.sort(reverse=True)
        best[c] = lst[0][1]
        if len(lst) > 1:
            straddle += 1

    # Tile origin from the minimum sampled coordinate.
    origins = {f: (min(tile_x[f]), min(tile_y[f])) for f in tile_x}

    # Row-major ordering over physical position.
    ys = sorted({round(o[1], 0) for o in origins.values()})
    row_of = {}
    for f, (ox, oy) in origins.items():
        row_of[f] = min(range(len(ys)), key=lambda i: abs(ys[i] - oy))
    ordered = sorted(origins, key=lambda f: (row_of[f], origins[f][0]))
    fov_id = {f: i + 1 for i, f in enumerate(ordered)}

    out = pd.DataFrame({
        "cell_id": list(best.keys()),
        "fov_name": [best[c] for c in best],
    })
    out["fov"] = out["fov_name"].map(fov_id).astype(int)
    out = out[["cell_id", "fov", "fov_name"]].sort_values("fov")

    meta = {
        "source": os.path.abspath(transcripts_path),
        "n_cells": int(len(out)),
        "n_tiles": int(len(fov_id)),
        "cells_straddling_tiles": int(straddle),
        "fraction_straddling": round(straddle / max(len(out), 1), 6),
        "ordering": ("row-major over measured tile origins; an assumed "
                     "acquisition order, no timestamp exists"),
        "fov_name_to_fov": {k: int(v) for k, v in sorted(fov_id.items(),
                                                         key=lambda kv: kv[1])},
    }
    return out, meta


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--transcripts", required=True,
                   help="vendor transcripts.parquet carrying fov_name")
    p.add_argument("--output", required=True,
                   help="output prefix; writes <prefix>.csv and <prefix>.json")
    p.add_argument("--quiet", action="store_true")
    a = p.parse_args()

    out, meta = build_map(a.transcripts, verbose=not a.quiet)
    out.to_csv(f"{a.output}.csv", index=False)
    with open(f"{a.output}.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"Wrote {a.output}.csv: {meta['n_cells']} cells over "
          f"{meta['n_tiles']} native tiles")
    print(f"  cells straddling more than one tile: "
          f"{meta['cells_straddling_tiles']} ({100 * meta['fraction_straddling']:.3f}%)")


if __name__ == "__main__":
    main()
