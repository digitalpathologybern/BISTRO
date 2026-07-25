# BISTRO

**B**ias **I**dentification in Image-based **S**patial **TR**anscript**O**mics

BISTRO is a Nextflow (DSL2) pipeline that benchmarks 11 normalization methods on
imaging-based spatial transcriptomics (iST) data — CosMx, Xenium, and MERSCOPE.
It runs every method in parallel, then evaluates each on FOV batch effects,
mean–variance behaviour, HVG selection stability, and cell-type-annotation
agreement, and finally bundles all metrics into a self-contained HTML report.

This repository accompanies the BISTRO manuscript and contains the full,
reproducible pipeline used to generate the results.

---

## Table of contents

1. [Pipeline overview](#pipeline-overview)
2. [Repository layout](#repository-layout)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Demo dataset](#demo-dataset)
6. [Input data](#input-data)
7. [Configuration](#configuration)
8. [Running the pipeline](#running-the-pipeline)
9. [Outputs](#outputs)
10. [Resuming / caching](#resuming--caching)
11. [License](#license)

---

## Pipeline overview

```
filtered.zarr ──► readZarr ──► assignFOV ──┬──► 11 normalizations (parallel)
                                            │
                                            ├──► HVG selection ──► UpSet plot
                                            │                  └─► pathway analysis
                                            │
                                            ├──► InSituType cell-type annotation
                                            │
                                            └──► Evaluation:
                                                  • FOV batch-effect (OLS / MELM)
                                                  • Transformation / variance stabilization
                                                  • HVG benchmark (per layer × fraction)
                                                  • Annotation agreement (ARI)
                                                          │
                                                          ▼
                                                  BISTRO_report.html
```

Normalizations evaluated: `none`, `cpm`, `cp10k`, `cp100`, `tmm`, `scran`,
`deseq2`, `areaNorm`, `SpaNorm-logpac`, `SpaNorm-pearson`, `scTransform`.

All spatial coordinates are in **micrometers (µm)** throughout the pipeline.

---

## Repository layout

```
.
├── run_BISTRO.nf              # Top-level Nextflow workflow
├── sample_configs/            # One .config per dataset
├── bin/
│   ├── preprocessing/         # FOV assignment (assign_fov.py)
│   ├── normalization/         # 11 normalization scripts (R + Python)
│   ├── annotation/            # InSituType (R)
│   ├── evaluation/            # Batch effect, transformation, HVG bench, ARI
│   ├── plotters/              # UpSet, pathway, GSVA, colormaps
│   ├── report/                # HTML report generator
│   ├── find_hvg.py
│   ├── save_expression_matrix.py
│   ├── reduce_dimensions.R
│   ├── do_clustering_banksy.R
│   └── do_clustering_lou_lei.R
├── utils/
│   ├── data_loader.py         # SpatialData I/O, layers, FOV merging, tissue annotations
│   └── helpers.py             # Stats helpers, FOV rasterization, layer-name parsing
├── envs/
│   ├── environment.yml        # Conda env, exact build pins
│   └── R_packages.csv         # All installed R packages and their versions
├── examples/
│   ├── run_bistro_local.sh    # Sample bash launcher (see below)
│   └── demo-dataset/          # Self-contained simulated demo (see below)
│       ├── make_demo_dataset.py
│       ├── make_demo_reference.R
│       ├── demo.config
│       ├── run_demo.sh
│       └── ground_truth.json
├── LICENSE
├── .gitignore
└── README.md
```

---

## Requirements

| Component | Version | Notes |
|---|---|---|
| Nextflow | ≥ 23.04 | DSL2 |
| Java     | ≥ 17    | required by Nextflow |
| Python   | 3.11    | conda env provided |
| R        | ≥ 4.4   | with Bioconductor 3.20 |

The exact dependency manifests used to produce the manuscript results are
committed to this repo under `envs/`:

| File | What it captures |
|---|---|
| `envs/environment.yml` | full Conda env (`conda env export --no-builds`) — every package and version |
| `envs/R_packages.csv`  | all R packages installed in the manuscript environment, with versions |

---

## Installation

```bash
# 1. Clone
git clone https://github.com/digitalpathologybern/BISTRO.git
cd BISTRO

# 2. Install Nextflow (or load a cluster module: `module load Nextflow`)
curl -s https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/

# 3. Recreate the Python conda environment.
conda env create -n bistro -f envs/environment.yml
conda activate bistro

# 4. Install the R / Bioconductor packages listed in envs/R_packages.csv.
#    BiocManager will resolve CRAN + Bioconductor packages; GitHub-only
#    packages (InSituType, SpatialPCA, SeuratWrappers, scGSVA) must be
#    installed manually via remotes::install_github().
Rscript -e '
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  pkgs <- read.csv("envs/R_packages.csv", stringsAsFactors = FALSE)$Package
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
'
```
---

## Demo dataset

A self-contained, fully simulated dataset lives in
[`examples/demo-dataset/`](examples/demo-dataset/). It runs the entire pipeline
on a single machine in minutes and needs no access to the manuscript data.

```bash
conda activate bistro
bash examples/demo-dataset/run_demo.sh
```

That regenerates the inputs, builds the InSituType reference, and runs
BISTRO end to end, leaving the report at
`examples/demo-dataset/demo_output/BISTRO_report.html`. A verified run takes
about 9 minutes on 4 cores and peaks at 2.8 GB.

On an HPC cluster, submit it rather than running it on a login node:

```bash
sbatch examples/demo-dataset/run_demo_slurm.sbatch
```

The simulation is not arbitrary: it plants a **known ground truth** — a per-FOV
library-size offset of known variance, a systematic decline in that offset
across the acquisition order, a designated set of marker genes, and five cell
types — and writes the realised values to `ground_truth.json`. A demo run can
therefore be checked for correctness, not merely for completing.
[`examples/demo-dataset/README.md`](examples/demo-dataset/README.md) lists what
to compare against what.

The demo doubles as the executable specification of the input schema described
in the next section.

---

## Input data

BISTRO starts from data that has **already been quality-controlled**. The
pipeline performs no cell or gene filtering of its own: whatever is in
`tables['filtered']` is what gets benchmarked. Producing that archive from a
vendor export is out of scope for this repository.

Every dataset needs three files (four with a GeoJSON annotation):

| File | Required | Purpose |
|---|---|---|
| `filtered.zarr` | yes | QC'd expression + cell metadata (schema below) |
| scRNA-seq reference | if `skip_annotation = false` | InSituType cell-type calling |
| tissue annotation | yes | fixed effect in the batch-effect model |
| H&E alignment matrix | only with a GeoJSON annotation | maps H&E pixels to image pixels |

`examples/demo-dataset/make_demo_dataset.py` is a working, runnable
implementation of everything below — when this document and that script
disagree, the script is correct.

### 1. `filtered.zarr` — schema

A [SpatialData](https://spatialdata.scverse.org) zarr archive. Only the table
is read; images, shapes and points are ignored (shape centroids are consulted
only as a last-resort fallback for coordinates).

**Structure**

| Element | Required | Notes |
|---|---|---|
| `tables['filtered']` | **yes** | the AnnData that the whole pipeline operates on |
| `tables['table']` | no | if present, re-indexed alongside `filtered` for MERSCOPE `EntityID` data |

**`X` — the expression matrix**

* **Raw integer counts**, cells × genes. Dense or `scipy.sparse`; both are handled.
* **Must not be normalized or log-transformed.** Every method starts from these
  counts, and TMM, scran and DESeq2 estimate size factors that are only
  meaningful on raw counts.
* Negative-control probes should already be **removed** from `X`; their
  per-cell total is carried in `obs` instead (see `total_counts_Negative`).

**`var` — genes**

* `var_names` holds gene symbols and must be unique.
* Symbols must overlap the scRNA-seq reference's row names, or annotation
  produces no usable genes and the run stops.
* The R stages normalize `-`, `_`, `:`, `/` and spaces to `.`
  (`clean_gene_names()` in `bin/annotation/run_insitutype.R`, and the same
  substitutions in `bin/plotters/hvg_upset_plot.py`). Symbols free of those
  characters are safest.

**`obs` — the index**

Cell IDs, unique and stable. They are the join key between the expression
matrix, the enriched metadata, the tissue annotation and every annotation
output, so they must survive a CSV round-trip unchanged. If an `EntityID`
column is present (the MERSCOPE convention) it is promoted to the index
automatically.

**`obs` — required columns**

All coordinates are in **micrometers**; there is no unit conversion anywhere in
the pipeline.

| Column | Type | Used by |
|---|---|---|
| `x_local_um`, `y_local_um` | float | pseudo-FOV rasterization; SpaNorm coordinates on Xenium/MERSCOPE |
| `x_global_um`, `y_global_um` | float | SpaNorm coordinates on CosMx; GeoJSON point-in-polygon; spatial plots |
| `area_um2` | float | `areaNorm` size factors — this method has no toggle, so the column is always required |

**`obs` — conditionally required**

| Column | When | Notes |
|---|---|---|
| `total_counts_Negative` | `skip_annotation = false` | **Sum** of negative-probe counts per cell. `run_insitutype.R` divides it by 20 to get the per-cell background mean. |
| `fov` | platforms with native FOVs (CosMx) | Integer FOV index. If absent, or constant, `assignFOV` rasterizes pseudo-FOVs from `x_local_um` / `y_local_um` using the technology's tile size. |
| `fov_center_x_um`, `fov_center_y_um` | never supplied by hand | Computed by `assignFOV` when missing. |

Any other columns are carried through to the metadata CSV untouched.

> **FOV numbering is the acquisition order.** The drift analysis regresses the
> per-FOV random intercepts against the FOV index, so for native-FOV platforms
> `fov` should increase in the order the instrument imaged the tiles.
> Rasterized platforms get a row-primary `fov` plus a column-primary `fov_perp`
> so the two scan directions can be compared.

### 2. scRNA-seq reference

Two accepted forms, selected by file extension:

| Extension | Expected contents |
|---|---|
| `.RData` | an object named exactly **`profile_matrix`** — genes × cell types, mean expression on the linear scale |
| `.rds` | a Seurat object with a `CellType` column; profiles are built with `AggregateExpression` |

### 3. Tissue annotation

Required: the batch-effect model's fixed effect is
`log_LS ~ C(tissue_annotations)`, so there is no path through
`bin/evaluation/batch_effect.py` without it.

| Format | Requirements |
|---|---|
| `.csv` | a `tissue_annotations` column (`niche` and `banksy_0.8` are accepted and renamed), keyed by a `cell_ID` column or by the unnamed index column, matching `obs` index values |
| `.geojson` | a QuPath export whose `classification.name` gives the region label; assigned by point-in-polygon against `x_global_um` / `y_global_um` |

### 4. H&E alignment matrix

Only used with a GeoJSON annotation: a 3×3 affine as a headerless CSV, applied
to the polygons before they are scaled from image pixels to micrometers by
`pixelSize`.

---

## Configuration

Each dataset is described by one `.config` file in `sample_configs/`.
Use one of those as a template. The required parameters are:

```groovy
params {
    // ── Input data ──────────────────────────────────────────────────────
    zarrFile          = '/path/to/dataset_filtered.zarr'
    scReferenceFile   = '/path/to/scRNAseq_reference.RData'
    tissueAnnotation  = '/path/to/annotations.geojson'    // or .csv, or ''
    heAlignmentPath   = '/path/to/he_imagealignment.csv'  // Xenium only, or ''
    pixelSize         = 0.2125                            // µm per pixel

    // ── Dataset metadata ────────────────────────────────────────────────
    datasetName  = 'MyDataset'
    technology   = 'Xenium'        // 'CosMx' | 'Xenium' | 'MERSCOPE'

    // ── Platform-specific flags ─────────────────────────────────────────
    separate_fovs = '1'            // SpaNorm per-FOV ('1') or global ('0')
    add_global    = '0'

    // ── Output ──────────────────────────────────────────────────────────
    output_folder_path = '/path/to/output/MyDataset/'
    outputNorm = 'norm'; outputHVG = 'hvg'; outputAnno = 'annotation'
    outputEval = 'evaluation'; outputPlots = 'plots'; outputPathway = 'pathway'

    // ── HVG benchmark ───────────────────────────────────────────────────
    hvg_fractions        = '0.1 0.25 0.5 0.75 1.0'
    hvg_n_bootstrap      = 100
    hvg_subsample_frac   = 0.80
    hvg_n_random_repeats = 5
    hvg_k_neighbors      = 15
    hvg_resolution       = 0.5
    hvg_seed             = 42
    hvgThreshold         = '0.75'

    n_boot_ci            = 200     // bootstrap reps for τ² CI

    // ── Reference annotation (HVG benchmark, Phase 4b) ───────────────────
    // Optional. Leave empty to use the pipeline's own InSituType output for
    // the `none` normalization at 100% HVG. Must be set to an existing
    // annotation CSV when skip_annotation = true and skip_hvg_bench = false.
    referenceAnnotation = ''

    // ── Toggles ─────────────────────────────────────────────────────────
    skip_annotation = false
    skip_pathway    = false
    skip_hvg_bench  = false
    restore_published = true       // skip a step if its outputs already exist
}
```

### Reference annotation for the HVG benchmark

Phase 4b of the HVG benchmark (marker-gene AUROC) scores each normalization
against a fixed set of cell-type labels. By default those labels are the
InSituType calls for the `none` normalization at 100% HVG
(`annotation/1.0_HVG_none_annotation.csv`), which the pipeline produces itself.

Because that file only exists when the annotation step runs, the combination
`skip_annotation = true` + `skip_hvg_bench = false` requires an explicit
reference:

```groovy
referenceAnnotation = '/path/to/1.0_HVG_none_annotation.csv'
```

The file must have cell IDs in the index column and labels in a `sup.clust`
column (a single-column CSV also works). Setting it when the annotation step
*does* run overrides the derived reference, which is useful for scoring
against an external ground truth. If the combination is misconfigured, the
run aborts before any process is scheduled with a message naming the three
ways out.

The `process { ... }` block at the bottom of the file controls resources and
`module load` commands; adjust for your cluster.

---

## Running the pipeline

### Local / single-machine

```bash
nextflow run run_BISTRO.nf -c sample_configs/xenium-cancerBreast-5k.config
```

### Reproducible bash launcher

A ready-to-edit launcher is in [`examples/run_bistro_local.sh`](examples/run_bistro_local.sh).
Copy it, edit the four paths at the top, and run:

```bash
cp examples/run_bistro_local.sh my_run.sh
# edit BISTRO_DIR, CONFIG_PATH, WORK_DIR, NXF_CACHE
bash my_run.sh
```

The launcher pins per-run `workDir` and Nextflow cache so multiple datasets do
not collide, and it forwards exit codes so it is safe to wrap in `cron`,
`tmux`, or your own scheduler.

---

## Outputs

All outputs are written to `params.output_folder_path` and mirror the
sub-directory names declared in the config:

```
<output_folder_path>/
├── norm/                                # one CSV per normalization (cells × genes)
│   ├── *_counts.csv, *_metadata.csv     # raw matrix + FOV-enriched metadata
│   ├── *_areaNorm.csv, *_cpm.csv, …     # 11 normalization layers
│   └── *_areaSF.csv, *_cpmSF.csv, …     # size factors (where applicable)
├── hvg/                                 # per-layer HVG tables + top-fraction selections
├── annotation/                          # InSituType cell-type calls per HVG fraction
├── pathway/                             # GO / MSigDB enrichments on HVGs
├── plots/                               # UpSet, pathway, GSVA, spatial plots
└── evaluation/
    ├── batch_effect/                    # OLS + MELM summaries, random intercepts, FOV summary
    ├── transformation/                  # mean–variance, PC1 correlation, transform summary
    ├── hvg_benchmark/                   # per-layer × per-fraction clustering / coherence / AUROC
    └── annotation_agreement/            # pairwise ARI between normalizations
└── BISTRO_report.html                   # self-contained final report
```

---

## Resuming / caching

The pipeline caches as:

**`restore_published = true`** (default) — each process checks
   `params.output_folder_path/<step>/` *before* running and short-circuits via
   symlinks if the expected output already exists. Set `restore_published =
   false` in the config to force a clean recomputation.

To rerun a single step from scratch, delete its sub-directory under the output
folder.

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

BISTRO orchestrates third-party normalization methods that carry their own
licenses (among them edgeR, DESeq2, scran, SpaNorm, Seurat/SCTransform and
InSituType). Each is invoked as a separate process rather than linked, but if
you redistribute BISTRO together with those dependencies, check their terms.

---
