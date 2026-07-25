# BISTRO demo dataset

A small, fully simulated imaging-based spatial transcriptomics dataset that
runs the complete BISTRO pipeline on one machine in minutes, with **no access
to the manuscript data required**.

It serves two purposes:

1. **A runnable demo**: proves the installation works end to end.
2. **A correctness check**: the simulation plants a known ground truth, so a
   run can be verified against what was injected rather than merely observed to
   finish.

It also acts as the executable specification of the `filtered.zarr` input
schema documented in the [top-level README](../../README.md#input-data): if the
prose and `make_demo_dataset.py` ever disagree, the script is correct.

---

## Running it

```bash
conda activate bistro          # environment from ../../envs/environment.yml
bash run_demo.sh
```

`run_demo.sh` does three things: generates the simulated inputs, converts the
reference profiles into the `.RData` InSituType expects, and runs the pipeline
with `demo.config`.

To regenerate only the data:

```bash
python make_demo_dataset.py --output_dir .
Rscript make_demo_reference.R demo_reference_profiles.csv demo_reference.RData
```

Generation is deterministic given `--seed` (default 42), so the same seed
always yields byte-identical inputs.

Results land in `demo_output/`, with the report at
`demo_output/BISTRO_report.html`.

### On an HPC cluster

`run_demo.sh` runs everything in the foreground, which is fine on a laptop but
should not be done on a login node. On the University of Bern cluster, submit
it as a single job instead:

```bash
sbatch examples/demo-dataset/run_demo_slurm.sbatch
```

That allocates one node, generates the data, runs the pipeline with
`demo.cluster.config` layered on top of `demo.config` (which adds the
`module load` blocks per process label), and finishes by verifying the
outputs. Adjust `--partition` and `R_LIBS_USER` for your site.

The overlay keeps `demo.config` itself free of any cluster assumptions, so the
same config works unchanged off-cluster.

### Verifying a run

```bash
python check_demo.py --demo_dir . --output_dir demo_output
```

This compares the run against `ground_truth.json` and exits non-zero if any
check fails, so it can be used as a regression test. It distinguishes two
kinds of result:

* `[CHECK]` is a principled pass/fail: a hypothesis test, the sign of an
  effect, or whether a confidence interval covers the injected value. These
  determine the exit code.
* `[INFO]` is a magnitude whose expected range has to be established from a
  reference run. Reported, never failed on.

The checks are described in the table below.

### Reference run

For calibration, a verified run of the committed inputs (seed 42, 4 cores,
16 GB, Nextflow 24.04.4, R 4.4.2):

| | |
|---|---|
| Wall time | 9 min 03 s |
| Peak memory | 2.8 GB |
| Nextflow tasks | 54 |
| Report | 1.1 MB, self-contained HTML |
| Verification | 13 checks passed, 0 failed |

Headline numbers from that run:

| Quantity | Injected | Recovered |
|---|---|---|
| FOV variance τ² | 0.04580 | 0.04084 (ratio 0.89) |
| Drift slope per FOV | −0.02331 | −0.02190 |
| Drift correlation | r = −0.519 | r = −0.497 |
| Cell types | 5 | ARI = 1.000 |
| Marker genes in top 25 % HVG | 40 planted | 26 recovered, 10 expected by chance (p = 1.2 × 10⁻⁹) |
| Mean–variance slope, raw counts | NB, θ = 10 | 1.288 (→ −0.031 under log(x+0.5)) |
| LRT for the FOV effect | n/a | p ≈ 0, ΔBIC = 1612.5 |

Exact values will shift slightly with BLAS threading and package versions;
the checks are written to tolerate that.

---

## What gets generated

| File | Contents |
|---|---|
| `demo_filtered.zarr` | SpatialData archive, table `filtered`, the pipeline input |
| `demo_tissue_annotations.csv` | tissue region per cell, the batch-effect model's fixed effect |
| `demo_reference_profiles.csv` | genes × cell types, mean expression |
| `demo_reference.RData` | the above as `profile_matrix`, for InSituType |
| `demo_color_list.txt` | colour list for `create_colormap.R` |
| `ground_truth.json` | simulation parameters **and realised values** |

The inputs are committed to the repository as well as being reproducible from
the generator, so the demo can be inspected without running anything.
`demo_reference.RData` is the one exception: it is rebuilt from
`demo_reference_profiles.csv` by `make_demo_reference.R`, since generating it
requires R.

`demo_output/`, the pipeline's results, is deliberately **not** committed. It
is fully reproducible from the inputs above, it is two orders of magnitude
larger than them, and it contains machine-specific paths, so a committed copy
would generate noisy diffs and risk drifting out of sync with the code.

---

## The simulation

2000 cells × 200 genes, on a 4 × 4 grid of 510.72 µm FOVs (the CosMx
convention, so `assignFOV` takes the native-FOV path), across 3 tissue regions
and 5 cell types.

Tissue regions are **concentric rings** about the slide centre (a tumour core,
a stromal rim, an immune periphery) rather than bands. This matters more than
it looks. The batch-effect model is
`log_LS ~ C(tissue_annotations) + (1|FOV)`, and the drift analysis regresses
the per-FOV random intercepts on the FOV index. Horizontal bands are nearly
collinear with a row-major FOV index (measured at r = −0.84 in an earlier
version of this simulation), so the tissue fixed effect absorbs the injected
acquisition drift and the drift analysis finds almost nothing. Radial distance
is symmetric about the middle rows, so it is nearly uncorrelated with
acquisition order (r ≈ 0.00) while staying spatially contiguous, which the
HVG spatial-coherence phase needs.

Counts are negative-binomial,

```
log s_i = u_f(i) + e_i                       per-cell log size factor
u_f     = DRIFT_SLOPE · centred_fov_index + N(0, FOV_NOISE_SD²)
mu_ig   = s_i · exp(beta[k(i), g] + gamma[t(i), g])
y_ig    ~ NB(mean = mu_ig, dispersion = theta)
```

which reproduces the `Var = mu + mu²/theta` mean–variance relationship that the
transformation analysis is built to characterise, at a realistic ~40 % zero
fraction.

Four things are deliberately planted:

| Injected | Recovered by |
|---|---|
| per-FOV offset `u_f`, of known variance | the MELM random-intercept variance τ² (Equation 2) |
| systematic decline in `u_f` across FOV index | `drift_slope` / `drift_pearson_r` |
| 40 cell-type marker genes carrying the between-type variance | HVG selection |
| 5 cell types with distinct profiles | InSituType annotation |

Gene symbols are real and free of `-`, `_`, `:`, `/` and spaces, so they pass
through the R stages' name normalization unchanged.

---

## Checking a run against the ground truth

`ground_truth.json` holds both the nominal parameters and the **realised**
values for the generated draw. The realised values are what a run should
recover.

| Check | Ground truth | Pipeline output |
|---|---|---|
| FOV batch effect is detected | `batch_effect.u_fov` is non-zero by construction | `lrt_pvalue` in `evaluation/batch_effect/BISTRO-Demo_batch_effect_summary.csv` should be far below 0.05, and `delta_bic` strongly positive, for the `none` layer |
| Magnitude of the batch effect | `batch_effect.realised_var_u_fov` | `var_random_intercept_reml`, same file, `none` layer, judged against the χ² sampling band for 16 FOVs (see below) |
| Acquisition drift | `batch_effect.realised_drift_slope` and `realised_drift_pearson_r` | `drift_slope` and `drift_pearson_r`, same file; sign must match and the correlation should track the realised one |
| Cell types | `obs['true_cell_type']` in the zarr | `annotation/1.0_HVG_none_annotation.csv`; labels are recovered only up to a permutation, so score with adjusted Rand index, not accuracy |
| Marker genes | `hvg.all_marker_genes` | the `none` column of `hvg/top_0.25_selected_hvg_per_method.csv` should be enriched for them relative to chance |
| Mean–variance behaviour | `nb_dispersion_theta` | the raw log–log slope in `evaluation/transformation/BISTRO-Demo_transformation_summary.csv` should sit near 1–2 (Poisson to NB) and fall towards 0 under variance-stabilising transformations |

**On τ².** The mixed-effects estimate is not expected to match
`realised_var_u_fov` exactly. It is a REML estimate from 16 groups, so it is
shrunk toward zero, and it is fitted on `log1p(library size)` rather than on
the latent log size factor. `check_demo.py` therefore tests it against the
sampling band that any variance estimate from *k* groups inherits,
`(k−1)/χ²_{0.975,k−1}` … `(k−1)/χ²_{0.025,k−1}`. For k = 16 that is a ratio of
0.55 to 2.40, rather than against an arbitrary tolerance.

**On drift: compare against the realised value, not the nominal one.** With
only 16 FOVs the drift actually present in a given draw scatters widely around
`DRIFT_SLOPE`. The generator therefore records `realised_drift_slope` and
`realised_drift_pearson_r`, the regression of the simulated offsets `u_fov` on
the FOV index, and those are what the pipeline can recover. The nominal
parameter is a property of the generative model; the realised value is a
property of the data on disk.

This also drove the choice of `DRIFT_SLOPE`. At an earlier value of −0.018 the
drift was comparable to the FOV noise, and the seed-42 draw landed at
r = −0.16, in the bottom 4 % of draws, leaving the drift analysis almost
nothing to find. The current value makes the systematic component dominant, so
the demonstration holds across seeds rather than depending on a lucky one.

**On the bootstrap CI.** `var_ci_lower` … `var_ci_upper` comes from
`bootstrap_var_ci()`, which resamples cells *within* each FOV while holding the
FOV set fixed. It is therefore a **conditional** interval: the uncertainty in
τ² given these 16 FOVs. That is the right question for characterising a single
slide, which is what BISTRO is for, but it is *not* a confidence interval for
the population variance of the FOV offsets, because it carries no uncertainty
about which 16 offsets were drawn. At k = 16 a population interval would span
roughly ±70 %, while this one spans a few percent, so it will not generally
cover the injected value. `check_demo.py` reports it without failing on it.

---

## How the demo config differs from the manuscript settings

`demo.config` is tuned for speed on a laptop and is **not** a template for real
runs. Copy from [`../../sample_configs/`](../../sample_configs/) instead.

| Parameter | Demo | Manuscript |
|---|---|---|
| `hvg_fractions` | `0.25 1.0` | `0.1 0.25 0.5 0.75 1.0` |
| `hvg_n_bootstrap` | 10 | 100 |
| `hvg_n_random_repeats` | 2 | 5 |
| `n_boot_ci` | 20 | 200 |
| `skip_pathway` | `true` | `false` |
| `restore_published` | `false` | `true` |
| executor | `local` | `slurm` |

`1.0` has to stay in `hvg_fractions`: the HVG benchmark's reference labels are
the InSituType calls at 100 % HVG on the `none` normalization, and
`hvgThreshold` must name one of the fractions actually computed.

Pathway analysis is off because `bin/plotters/pathway_analysis_HVG.R` needs
`GOfuncR` and `msigdbr`, and the `pathway_analysis_hvg` process issues a
`module load CMake` that only exists on an HPC with environment modules.

SpaNorm runs in global rather than per-FOV mode (`separate_fovs = '0'`): at
~125 cells per FOV the demo sits below the 300-cell threshold at which
`spanorm.R` merges FOVs for spline stability, so per-FOV mode would collapse
them into a single group anyway.

---

## Adapting it

The constants at the top of `make_demo_dataset.py` control the simulation:
`N_CELLS`, `N_GENES`, `FOV_GRID`, `DRIFT_SLOPE`, `FOV_NOISE_SD`, `THETA`,
`MARKER_LOG_FC` and the cell-type and tissue-region definitions.

Two are worth understanding before changing them:

* **`FOV_GRID`** sets the number of groups the mixed-effects model estimates
  τ² from. Fewer FOVs means a noisier, more heavily shrunk estimate.
* **`THETA`** is the negative-binomial dispersion. Lowering it increases
  overdispersion and makes variance stabilisation harder, which is useful for
  stress-testing the transformation analysis.

To simulate a rasterized platform instead, set `technology` to `Xenium` or
`MERSCOPE` in `demo.config` and drop the `fov` column from `obs`; `assignFOV`
will then generate both row-primary and column-primary pseudo-FOVs and the
report will include the scan-direction drift comparison.
