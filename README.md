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
5. [Input data](#input-data)
6. [Configuration](#configuration)
7. [Running the pipeline](#running-the-pipeline)
8. [Outputs](#outputs)
9. [Resuming / caching](#resuming--caching)
10. [Citation](#citation)

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
│   └── run_bistro_local.sh    # Sample bash launcher (see below)
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

## Input data

For every dataset you need:

1. **`filtered.zarr`** — A [SpatialData](https://spatialdata.scverse.org) zarr
   archive with a quality-controlled `tables['filtered']` AnnData. Cells must
   have spatial coordinates in micrometers (column convention:
   `x_local_um`, `y_local_um`, `x_global_um`, `y_global_um`, `area_um2`).
2. **scRNA-seq reference** — `.rds` or `.RData` Seurat/SCE object used by
   InSituType for cell-type annotation.
3. **Tissue annotation (optional)** — `.csv` (with a `tissue_annotations`,
   `niche`, or `banksy_0.8` column) **or** `.geojson` (QuPath export). The
   pipeline will perform point-in-polygon assignment.
4. **H&E alignment matrix (optional)** — A 3×3 affine CSV (no
   header) mapping H&E pixels to pixels.

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
