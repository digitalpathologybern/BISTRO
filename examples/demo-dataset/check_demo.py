#!/usr/bin/env python
"""
BISTRO demo verification
========================

Checks a completed demo run against the ground truth planted by
make_demo_dataset.py, and reports whether BISTRO recovered what was injected.

Two kinds of check are reported:

  [CHECK]  A principled pass/fail criterion: a hypothesis test, a sign, or a
           confidence-interval coverage statement. These decide the exit code.

  [INFO]   A magnitude whose expected range has to be established empirically
           from a reference run. Reported, never failed on.

Exit code is 0 if every [CHECK] passes, 1 otherwise.

Usage
-----
    python check_demo.py --demo_dir . --output_dir demo_output
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2, hypergeom
from sklearn.metrics import adjusted_rand_score

DATASET = "BISTRO-Demo"
REF_LAYER = "none"          # the layer compared against ground truth
HVG_CHECK_FRACTION = "0.25"

# Layers every run should produce.
EXPECTED_LAYERS = {
    "none", "cpm", "cp10k", "cp100", "tmm", "scran", "deseq2", "areaNorm",
    "spanorm-logpac", "spanorm-pearson", "scTransform",
}


class Report:
    """Collects check outcomes and renders a summary."""

    def __init__(self):
        self.checks = []
        self.infos = []
        self.skips = []

    def check(self, name, passed, detail):
        self.checks.append((name, bool(passed), detail))
        flag = "PASS" if passed else "FAIL"
        print(f"  [CHECK] {flag}  {name}\n            {detail}")

    def info(self, name, detail):
        self.infos.append((name, detail))
        print(f"  [INFO]        {name}\n            {detail}")

    def skip(self, name, reason):
        self.skips.append((name, reason))
        print(f"  [SKIP]        {name}\n            {reason}")

    def summary(self):
        n_pass = sum(1 for _, p, _ in self.checks if p)
        n_fail = len(self.checks) - n_pass
        print()
        print("=" * 72)
        print(f"  {n_pass} passed, {n_fail} failed, "
              f"{len(self.infos)} informational, {len(self.skips)} skipped")
        if n_fail:
            print()
            for name, passed, detail in self.checks:
                if not passed:
                    print(f"  FAILED: {name} -- {detail}")
        print("=" * 72)
        return n_fail == 0


def _find(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


# ============================================================================
# Individual checks
# ============================================================================

def check_outputs_present(out_dir, rep):
    """Every stage should have produced its headline artifact."""
    expected = {
        "normalization matrices": f"{out_dir}/norm/*_counts_*.csv",
        "HVG tables": f"{out_dir}/hvg/*_hvg.csv",
        "HVG selections": f"{out_dir}/hvg/top_*_selected_hvg_per_method.csv",
        "annotations": f"{out_dir}/annotation/*_annotation.csv",
        "batch effect summary":
            f"{out_dir}/evaluation/batch_effect/*_batch_effect_summary.csv",
        "transformation summary":
            f"{out_dir}/evaluation/transformation/*_transformation_summary.csv",
        "HVG benchmark": f"{out_dir}/evaluation/hvg_benchmark/checkpoints/*.csv",
        "annotation agreement":
            f"{out_dir}/evaluation/annotation_agreement/*_pairwise_ari_*.csv",
        "HTML report": f"{out_dir}/BISTRO_report.html",
    }
    missing = [name for name, pat in expected.items() if not glob.glob(pat)]
    rep.check(
        "All pipeline stages produced output",
        not missing,
        "all present" if not missing else f"missing: {', '.join(missing)}",
    )


def check_layers(out_dir, rep):
    """All 11 normalizations should appear in the norm directory."""
    found = set()
    for path in glob.glob(f"{out_dir}/norm/*.csv"):
        base = os.path.basename(path)
        if any(tok in base for tok in ("metadata", "SF", "counts.csv")):
            continue
        found.add(os.path.splitext(base)[0].split("_")[-1])

    missing = EXPECTED_LAYERS - found
    rep.check(
        "All 11 normalization layers produced",
        not missing,
        f"found {len(found)}/{len(EXPECTED_LAYERS)}"
        + (f"; missing {sorted(missing)}" if missing else ""),
    )


def check_batch_effect(out_dir, gt, rep):
    """The injected per-FOV offset should be detected, signed correctly, and
    of the right magnitude."""
    path = _find(f"{out_dir}/evaluation/batch_effect/*_batch_effect_summary.csv")
    if path is None:
        rep.skip("Batch effect", "summary CSV not found")
        return

    df = pd.read_csv(path)
    row = df[df["layer"] == REF_LAYER]
    if row.empty:
        rep.skip("Batch effect", f"no '{REF_LAYER}' layer in {path}")
        return
    row = row.iloc[0]

    injected = gt["batch_effect"]["realised_var_u_fov"]

    # -- Detection: the FOV random intercept must be strongly supported -------
    p = float(row["lrt_pvalue"])
    rep.check(
        "FOV batch effect is detected (LRT)",
        p < 1e-3,
        f"lrt_pvalue = {p:.3e} (injected var(u_fov) = {injected:.5f})",
    )

    rep.check(
        "BIC prefers the FOV model",
        float(row["delta_bic"]) > 0,
        f"delta_bic = {float(row['delta_bic']):.1f}",
    )

    # -- Magnitude: tau^2 within the sampling band implied by k FOVs ----------
    # The estimand is the variance of the per-FOV offsets, estimated from
    # k = n_fovs groups. Even a perfect estimator inherits the sampling
    # variability of a variance estimate from k draws, whose 95% band is
    # (k-1)/chi2_{0.975,k-1} .. (k-1)/chi2_{0.025,k-1} times the truth. That
    # band -- not an arbitrary tolerance -- is the criterion.
    tau2 = float(row["var_random_intercept_reml"])
    k = int(gt["n_fovs"])
    lo_ratio = (k - 1) / chi2.ppf(0.975, k - 1)
    hi_ratio = (k - 1) / chi2.ppf(0.025, k - 1)
    ratio = tau2 / injected

    rep.check(
        "tau^2 is within the sampling band implied by the FOV count",
        lo_ratio <= ratio <= hi_ratio,
        f"REML {tau2:.5f} vs injected {injected:.5f} (ratio {ratio:.2f}); "
        f"expected band for k={k} FOVs is [{lo_ratio:.2f}, {hi_ratio:.2f}]",
    )

    # -- Bootstrap CI: reported, not failed on -------------------------------
    # bootstrap_var_ci() now defaults to scheme='cluster': FOVs are resampled
    # with replacement, so the interval DOES carry the uncertainty about which
    # k offsets were drawn and is the interval for tau^2 as a variance
    # component. That is what a cross-slide comparison needs, and it is what
    # the manuscript's claims are about. With k=16 it is wide, on the order of
    # +/-70%, and it is expected to cover the injected value most of the time.
    # scheme='within' (cells resampled inside a fixed FOV set) remains
    # available and gives the much narrower CONDITIONAL interval.
    lo, hi = float(row["var_ci_lower"]), float(row["var_ci_upper"])
    if np.isfinite(lo) and np.isfinite(hi):
        rep.info(
            "Conditional bootstrap CI for tau^2",
            f"[{lo:.5f}, {hi:.5f}] (FOV cluster resampling); "
            f"injected {injected:.5f} "
            f"{'inside' if lo <= injected <= hi else 'OUTSIDE'} the interval",
        )

    # -- Drift ---------------------------------------------------------------
    # Compared against the drift REALISED in this draw, not the nominal
    # DRIFT_SLOPE: with 16 FOVs the realised slope scatters appreciably around
    # the nominal one, and the realised value is what is actually present in
    # the data for the pipeline to find.
    be = gt["batch_effect"]
    nominal_slope = be["drift_slope_per_fov"]
    realised_slope = be.get("realised_drift_slope", nominal_slope)
    realised_r = be.get("realised_drift_pearson_r")
    got_slope = float(row["drift_slope"])
    got_r = float(row["drift_pearson_r"])

    rep.check(
        "Acquisition drift recovered with the injected sign",
        np.sign(got_slope) == np.sign(realised_slope),
        f"realised slope {realised_slope:+.5f}/FOV "
        f"(nominal {nominal_slope:+.4f}), recovered {got_slope:+.5f}",
    )

    if realised_r is not None:
        rep.check(
            "Recovered drift correlation tracks the realised one",
            abs(got_r - realised_r) < 0.25,
            f"realised r = {realised_r:+.3f}, recovered r = {got_r:+.3f} "
            f"(difference {abs(got_r - realised_r):.3f}; the random intercepts "
            f"are shrunk slightly and share variance with the tissue term, so "
            f"exact equality is not expected)",
        )


def check_annotation(out_dir, rep):
    """InSituType should recover the simulated cell types up to a permutation."""
    anno = _find(f"{out_dir}/annotation/1.0_HVG_{REF_LAYER}_annotation.csv")
    meta = _find(f"{out_dir}/norm/*_metadata.csv")
    if anno is None or meta is None:
        rep.skip("Cell-type annotation", "annotation or metadata CSV not found")
        return

    calls = pd.read_csv(anno, index_col=0)
    truth = pd.read_csv(meta, index_col=0)

    if "true_cell_type" not in truth.columns:
        rep.skip("Cell-type annotation", "true_cell_type absent from metadata")
        return

    col = "sup.clust" if "sup.clust" in calls.columns else calls.columns[0]
    calls.index = calls.index.astype(str)
    truth.index = truth.index.astype(str)
    shared = calls.index.intersection(truth.index)
    if len(shared) < 100:
        rep.skip("Cell-type annotation", f"only {len(shared)} shared cells")
        return

    ari = adjusted_rand_score(truth.loc[shared, "true_cell_type"].astype(str),
                             calls.loc[shared, col].astype(str))
    # ARI is 0 in expectation for independent labellings, so anything clearly
    # above 0 demonstrates real recovery. The threshold is deliberately loose;
    # tighten it once a reference run establishes the achievable value.
    rep.check(
        "Annotation recovers the simulated cell types",
        ari > 0.20,
        f"ARI = {ari:.3f} over {len(shared)} cells "
        f"(0 = chance, 1 = perfect)",
    )


def check_hvg_recovery(out_dir, gt, rep):
    """HVG selection should be enriched for the planted marker genes."""
    path = _find(
        f"{out_dir}/hvg/top_{HVG_CHECK_FRACTION}_selected_hvg_per_method.csv"
    )
    if path is None:
        rep.skip("HVG marker recovery", "HVG selection CSV not found")
        return

    sel = pd.read_csv(path, index_col=0)
    if REF_LAYER not in sel.columns:
        rep.skip("HVG marker recovery",
                 f"no '{REF_LAYER}' column in {os.path.basename(path)}")
        return

    # The pipeline normalizes -, _, :, / and spaces to '.' in gene names.
    def norm(g):
        for ch in "-_:/ ":
            g = g.replace(ch, ".")
        return g

    markers = {norm(g) for g in gt["hvg"]["all_marker_genes"]}
    universe = {norm(g) for g in sel.index}
    chosen = {norm(g) for g in sel.index[sel[REF_LAYER].astype(bool)]}

    markers &= universe
    n_total, n_markers = len(universe), len(markers)
    n_chosen, n_hit = len(chosen), len(chosen & markers)

    # One-sided hypergeometric: P(at least n_hit markers | random selection)
    p = hypergeom.sf(n_hit - 1, n_total, n_markers, n_chosen)
    expected = n_chosen * n_markers / n_total

    rep.check(
        "HVG selection is enriched for the planted markers",
        p < 0.01,
        f"{n_hit}/{n_markers} markers in the top {HVG_CHECK_FRACTION} "
        f"({n_chosen} genes); {expected:.1f} expected by chance; "
        f"hypergeometric p = {p:.3e}",
    )


def check_transformation(out_dir, rep):
    """Untransformed counts should show mean-variance dependence, and
    variance-stabilising transforms should reduce it."""
    path = _find(
        f"{out_dir}/evaluation/transformation/*_transformation_summary.csv"
    )
    if path is None:
        rep.skip("Transformation analysis", "summary CSV not found")
        return

    df = pd.read_csv(path)
    sub = df[df["normalization"] == REF_LAYER]
    if sub.empty:
        rep.skip("Transformation analysis", f"no '{REF_LAYER}' rows")
        return

    raw = sub[sub["transformation"] == "raw"]
    if raw.empty:
        rep.skip("Transformation analysis", "no untransformed row")
        return

    raw_slope = float(raw.iloc[0]["slope"])
    rep.check(
        "Raw counts show mean-variance dependence",
        raw_slope > 0.5,
        f"log-log slope = {raw_slope:.3f} "
        f"(0 = stabilized, 1 = Poisson, 2 = negative binomial)",
    )

    logged = sub[sub["transformation"].str.startswith("log(", na=False)]
    if not logged.empty:
        best = logged.loc[logged["slope"].abs().idxmin()]
        rep.check(
            "A log transform reduces mean-variance dependence",
            abs(float(best["slope"])) < abs(raw_slope),
            f"best is {best['transformation']} at slope "
            f"{float(best['slope']):.3f}, versus {raw_slope:.3f} raw",
        )


def check_annotation_agreement(out_dir, rep):
    """The pairwise ARI matrix should be square, symmetric and unit-diagonal."""
    path = _find(
        f"{out_dir}/evaluation/annotation_agreement/*_pairwise_ari_*.csv"
    )
    if path is None:
        rep.skip("Annotation agreement", "pairwise ARI CSV not found")
        return

    ari = pd.read_csv(path, index_col=0)
    square = ari.shape[0] == ari.shape[1]
    diag_ok = square and np.allclose(np.diag(ari.values), 1.0)
    vals = ari.values[~np.eye(ari.shape[0], dtype=bool)] if square else []
    symmetric = square and np.allclose(ari.values, ari.values.T,
                                       equal_nan=True)

    rep.check(
        "Pairwise ARI matrix is well formed",
        square and diag_ok and symmetric,
        f"{ari.shape[0]}x{ari.shape[1]}, unit diagonal: {diag_ok}, "
        f"symmetric: {symmetric}",
    )
    if len(vals):
        rep.info(
            "Between-method annotation agreement",
            f"off-diagonal ARI: median {np.nanmedian(vals):.3f}, "
            f"range [{np.nanmin(vals):.3f}, {np.nanmax(vals):.3f}]",
        )


def check_no_degenerate_outputs(out_dir, rep):
    """Guard against silently empty or all-NaN evaluation tables."""
    bad = []
    for path in glob.glob(f"{out_dir}/evaluation/**/*.csv", recursive=True):
        try:
            df = pd.read_csv(path)
        except Exception:
            bad.append(f"{os.path.basename(path)} (unreadable)")
            continue
        if len(df) == 0:
            bad.append(f"{os.path.basename(path)} (empty)")
            continue
        numeric = df.select_dtypes(include=[np.number])
        if len(numeric.columns) and numeric.isna().all().all():
            bad.append(f"{os.path.basename(path)} (all NaN)")

    rep.check(
        "No empty or all-NaN evaluation tables",
        not bad,
        "all tables populated" if not bad else "; ".join(bad),
    )


# ============================================================================
# Entry point
# ============================================================================

def main():
    p = argparse.ArgumentParser(
        description="Verify a BISTRO demo run against the injected ground truth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--demo_dir", default=os.path.dirname(os.path.abspath(__file__)),
                   help="Directory holding ground_truth.json")
    p.add_argument("--output_dir", default=None,
                   help="Pipeline output directory "
                        "(default: <demo_dir>/demo_output)")
    args = p.parse_args()

    out_dir = args.output_dir or os.path.join(args.demo_dir, "demo_output")
    gt_path = os.path.join(args.demo_dir, "ground_truth.json")

    if not os.path.exists(gt_path):
        print(f"ERROR: ground truth not found at {gt_path}", file=sys.stderr)
        print("Run make_demo_dataset.py first.", file=sys.stderr)
        sys.exit(2)
    if not os.path.isdir(out_dir):
        print(f"ERROR: output directory not found at {out_dir}", file=sys.stderr)
        print("Run the pipeline first (see run_demo.sh).", file=sys.stderr)
        sys.exit(2)

    with open(gt_path) as fh:
        gt = json.load(fh)

    print("=" * 72)
    print(f"  BISTRO demo verification")
    print(f"  ground truth : {gt_path}")
    print(f"  run output   : {out_dir}")
    print(f"  simulation   : {gt['n_cells']} cells x {gt['n_genes']} genes, "
          f"{gt['n_fovs']} FOVs, seed {gt['seed']}")
    print("=" * 72)
    print()

    rep = Report()

    print("-- Structural --------------------------------------------------")
    check_outputs_present(out_dir, rep)
    check_layers(out_dir, rep)
    check_no_degenerate_outputs(out_dir, rep)

    print()
    print("-- Batch effect (Section 2.2) ----------------------------------")
    check_batch_effect(out_dir, gt, rep)

    print()
    print("-- Transformation (Section 2.3) --------------------------------")
    check_transformation(out_dir, rep)

    print()
    print("-- HVG selection (Section 2.4) ---------------------------------")
    check_hvg_recovery(out_dir, gt, rep)

    print()
    print("-- Annotation (Section 2.4) ------------------------------------")
    check_annotation(out_dir, rep)
    check_annotation_agreement(out_dir, rep)

    ok = rep.summary()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
