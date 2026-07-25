# BISTRO

**B**ias **I**dentification in Image-based **S**patial **TR**anscript**O**mics

BISTRO is a Nextflow (DSL2) pipeline that benchmarks 11 normalization methods on
imaging-based spatial transcriptomics (iST) data: CosMx, Xenium, and MERSCOPE.
It runs every method in parallel, then evaluates each on FOV batch effects,
mean–variance behaviour, HVG selection stability, and cell-type-annotation
agreement, and finally bundles all metrics into a self-contained HTML report.

This repository accompanies the BISTRO manuscript and contains the full,
reproducible pipeline used to generate the results.

---

## Table of contents

1. [Pipeline overview](#pipeline-overview)
2. [Repository layout](#repository-layout)
3. [System requirements](#system-requirements)
4. [Installation](#installation)
5. [Demo](#demo)
6. [Input data](#input-data)
7. [Configuration](#configuration)
8. [Instructions for use](#instructions-for-use)
9. [Outputs](#outputs)
10. [Resuming / caching](#resuming--caching)
11. [Reproducing the manuscript results](#reproducing-the-manuscript-results)
12. [License](#license)

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
│   ├── requirements-lock.txt         # Exact pins, the install path
│   ├── environment.yml               # Readable list of direct dependencies
│   ├── install_R_packages.R          # R dependency installer and checker
│   └── R_packages.csv                # All 371 R packages, provenance only
├── examples/
│   ├── run_bistro_local.sh    # Sample bash launcher (see below)
│   └── demo-dataset/          # Self-contained simulated demo (see below)
│       ├── make_demo_dataset.py
│       ├── make_demo_reference.R
│       ├── check_demo.py
│       ├── demo.config
│       ├── demo.cluster.config
│       ├── run_demo.sh
│       ├── run_demo_slurm.sbatch
│       └── ground_truth.json
├── LICENSE
├── .gitignore
└── README.md
```

---

## System requirements

### Operating system

BISTRO is a Nextflow pipeline and has no OS-specific code, but it has only
been tested on Linux.

| | |
|---|---|
| **Tested on** | Rocky Linux 9.7 (Blue Onyx), kernel 5.14.0, x86_64 |
| Expected to work | Any x86_64 Linux with the dependencies below; macOS (Intel or Apple silicon), untested |
| Windows | Not supported directly. Use WSL2 with a Linux distribution, untested |

### Core software

| Component | Minimum | Tested with |
|---|---|---|
| Nextflow | 23.04 (DSL2) | **24.04.4** |
| Java | 11 | **11.0.20** and **17.0.6** |
| Python | 3.11 | **3.11.8** |
| R | 4.4 | **4.4.2** |
| Bioconductor | 3.20 | **3.20** |

### Python packages

Installed from [`envs/requirements-lock.txt`](envs/requirements-lock.txt),
which pins all 113 packages including transitive dependencies.
[`envs/environment.yml`](envs/environment.yml) is a readable summary of the
~20 packages BISTRO imports directly.

| Package | Version | | Package | Version |
|---|---|---|---|---|
| numpy | 1.26.4 | | scanpy | 1.10.2 |
| pandas | 2.0.0 | | statsmodels | 0.14.2 |
| scipy | 1.13.1 | | scikit-learn | 1.5.0 |
| anndata | 0.10.8 | | matplotlib | 3.9.0 |
| spatialdata | 0.2.6 | | geopandas | 1.0.0 |
| zarr | 2.15.0 | | shapely | 2.0.4 |
| numcodecs | 0.12.1 | | python-igraph | 0.11.8 |
| xarray | 2024.11.0 | | leidenalg | 0.10.2 |
| tqdm | 4.66.4 | | upsetplot | 0.9.0 |

> **On pandas, and why the install needs `--no-deps`.** The environment that
> produced the manuscript results holds pandas 2.0.0 together with
> xarray 2024.11.0, and xarray 2024.11.0 declares `pandas>=2.1`. It runs
> correctly, because the spatialdata code paths BISTRO exercises never touch
> the pandas 2.1 API, but pip will not construct that combination: every
> xarray at or above the `>=2024.10.0` floor spatialdata requires wants
> pandas 2.1 or newer. Since pandas 2.0.0 is the version every published
> result was computed on, it is held and pip's resolver is bypassed with
> `--no-deps`. This is why `requirements-lock.txt` must list every transitive
> dependency, and why it should be regenerated with `pip freeze` rather than
> edited by hand.

### R packages

Full manifest of all 371 packages in
[`envs/R_packages.csv`](envs/R_packages.csv). The ones that matter:

| Package | Version | Role |
|---|---|---|
| SpaNorm | 1.0.0 | spatially aware normalization |
| InSituType | 2.0 | cell-type annotation |
| scran | 1.34.0 | pooling-based size factors |
| DESeq2 | 1.46.0 | median-of-ratios size factors |
| edgeR | 4.4.2 | TMM size factors |
| Seurat | 5.5.0 | SCTransform |
| sctransform | 0.4.3 | regularized negative binomial |
| SingleCellExperiment | 1.28.1 | data structure |
| SpatialExperiment | 1.16.0 | data structure |
| limma | 3.62.2 | linear modelling |
| GOfuncR | 1.26.0 | pathway analysis (optional) |
| msigdbr | 26.1.0 | gene sets (optional) |

### System libraries

Several R packages link against C/C++ libraries that must be present at
**load** time, not just at install time. When one is missing the R package
installs fine and then fails to load, which is easy to misread as a missing
package.

| Library | Needed by | Debian/Ubuntu | Symptom if absent |
|---|---|---|---|
| GLPK | igraph → Seurat, scran | `libglpk-dev` | `libglpk.so.40: cannot open shared object file` |
| ImageMagick (Magick++) | magick → SpatialExperiment | `libmagick++-dev` | `libMagick++-7...: cannot open shared object file` |
| libxml2 | XML, xml2 | `libxml2-dev` | fails at install |
| GDAL, PROJ, GEOS | sf, terra | `libgdal-dev libproj-dev libgeos-dev` | fails at install |
| UDUNITS | units → sf | `libudunits2-dev` | fails at install |
| ICU | stringi | `libicu-dev` | fails at install |

Check all R dependencies, and distinguish "not installed" from "installed but
not loadable", with:

```bash
Rscript envs/install_R_packages.R --check
```

On a module-based HPC these are usually just modules you have not loaded; see
[`examples/demo-dataset/demo.cluster.config`](examples/demo-dataset/demo.cluster.config)
for the exact set used here.

### Hardware

**No non-standard hardware is required.** There is no GPU code path, and no
dependency on any accelerator, interconnect, or specialised storage. Any
x86_64 machine will do.

Memory is the binding constraint, and it scales with cells × genes:

| Workload | Cores | Memory | Disk |
|---|---|---|---|
| **Demo** (2,000 cells × 200 genes) | 2–4 | **8 GB** is ample; measured peak 2.8 GB | ~200 MB |
| Small panel (~10⁵ cells × 300–1,000 genes) | 4–8 | 32–64 GB | ~50 GB |
| Manuscript datasets (up to ~10⁶ cells × 18,000 genes) | 6+ | **256 GB**, and **768 GB** for scran and SpaNorm | 100–400 GB per dataset |

The demo runs comfortably on any modern laptop. The manuscript datasets do
not: they were run on an HPC cluster, and the memory figures above are the
values requested in [`sample_configs/`](sample_configs/). Normalized matrices
are written as dense CSV, which is what drives the disk figures. A single
18k-plex layer is 20–32 GB.

---

## Installation

There is nothing to compile: BISTRO is a Nextflow workflow over Python and R
scripts. Installation is entirely a matter of providing the interpreters and
their packages.

### 1. Clone

```bash
git clone https://github.com/digitalpathologybern/BISTRO.git
cd BISTRO
```

### 2. Nextflow and Java

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/
```

On a cluster, `module load Nextflow` instead. Nextflow needs Java 11 or newer
already on `PATH`.

### 3. System libraries

Install the libraries listed under
[System libraries](#system-libraries) *before* the R packages, since several will
not build without them. On Debian/Ubuntu:

```bash
sudo apt-get install -y libglpk-dev libmagick++-dev libxml2-dev \
    libgdal-dev libproj-dev libgeos-dev libudunits2-dev libicu-dev
```

### 4. Python environment

```bash
conda create -n bistro -c conda-forge python=3.11.8 setuptools=69.1.0 numpy=1.26.4
conda activate bistro
pip install --no-deps -r envs/requirements-lock.txt
```

**`--no-deps` is required, not optional.** It bypasses pip's dependency
resolver, which is the only way to reproduce the environment the manuscript
results were computed on; see the note on pandas under
[Python packages](#python-packages). Because `--no-deps` installs exactly the
lines in the lock file and nothing else, that file lists all 113 packages
including transitive ones.

### 5. R packages

```bash
Rscript envs/install_R_packages.R
```

This installs the pipeline's direct dependencies from CRAN and Bioconductor
and pulls InSituType from GitHub; CRAN and Bioconductor resolve the transitive
dependencies. Add `--optional` to include `scGSVA`, needed only for the GSVA
step that is disabled by default.

`envs/R_packages.csv` records all 371 packages present in the environment that
produced the manuscript results. It is a provenance record, not an install
list. Installing the direct dependencies is faster and less brittle.

### 6. Verify

```bash
Rscript envs/install_R_packages.R --check
```

This distinguishes "not installed" from "installed but will not load", the
latter being a missing system library rather than a missing R package. Then
run the [demo](#demo), which exercises every stage end to end.

### Typical install time

On a normal desktop with a broadband connection:

| Step | Time |
|---|---|
| Clone | seconds |
| Nextflow + Java | 1-2 min |
| System libraries | 1-3 min |
| Python environment (conda + pip) | **2 min 06 s**, measured |
| R packages | **1-3 hours** |
| **Total** | **roughly 1.5-3.5 hours**, nearly all of it R |

The Python step is measured: 2 min 06 s on 4 cores, producing a 1.5 GB
environment, verified by importing every module the pipeline uses, generating
the demo dataset, and reading it back through the pipeline's own loader.

The R step dominates and varies enormously with what CRAN and Bioconductor
ship as binaries for your platform. On Linux, `BiocManager` builds most
packages from source, and Seurat, DESeq2 and their dependency trees are the slow
part. On macOS and Windows, binaries are usually available and the step drops
to 15–30 minutes. If you already have a Bioconductor 3.20 installation, it can
be minutes.

---

## Demo

A self-contained, fully simulated dataset lives in
[`examples/demo-dataset/`](examples/demo-dataset/). It exercises every stage of
the pipeline on one machine in minutes and needs no access to the manuscript
data, which is far too large to serve as a demo.

The inputs are committed to the repository (~830 KB), so they can be inspected
without running anything, and they are also reproducible from the generator at
a fixed seed.

### Instructions

```bash
conda activate bistro
bash examples/demo-dataset/run_demo.sh
```

That runs four steps: generate the simulated inputs, build the InSituType
reference, run BISTRO, and verify the output against the injected ground
truth. It exits non-zero if any check fails.

On an HPC cluster, submit it as a job rather than running it on a login node:

```bash
sbatch examples/demo-dataset/run_demo_slurm.sbatch
```

That variant layers
[`demo.cluster.config`](examples/demo-dataset/demo.cluster.config) on top of
`demo.config` to load environment modules per process, and runs the whole demo
inside a single allocation.

### Expected run time

| | |
|---|---|
| **Wall time** | **9 min 03 s** |
| Hardware | 4 cores of an AMD EPYC 7742, 16 GB allocated |
| Peak memory | 2.8 GB |
| Nextflow tasks | 54 |
| Disk written | ~163 MB |

Measured on Rocky Linux 9.7 with the versions listed above. Peak memory stayed
well under the 16 GB allocated, so 8 GB is sufficient; expect a broadly
similar wall time on any current 4-core desktop.

Both entry points have been verified end to end on this dataset:
`run_demo.sh` and `run_demo_slurm.sbatch` each complete with 13 of 13 checks
passing.

### Expected output

Results land in `examples/demo-dataset/demo_output/`:

| Path | Contents |
|---|---|
| `BISTRO_report.html` | self-contained HTML report, ~1.1 MB |
| `norm/` | the 11 normalized count matrices |
| `hvg/` | per-method HVG tables and selections |
| `annotation/` | InSituType calls per normalization |
| `evaluation/batch_effect/` | τ², LRT, drift statistics |
| `evaluation/transformation/` | mean–variance slopes |
| `evaluation/hvg_benchmark/` | clustering stability, coherence, AUROC |
| `evaluation/annotation_agreement/` | pairwise ARI between normalizations |

The verification step prints a pass/fail line per check and ends with a
summary. A correct run reports **13 passed, 0 failed**:

```
  [CHECK] PASS  All 11 normalization layers produced
  [CHECK] PASS  FOV batch effect is detected (LRT)
            lrt_pvalue = 0.000e+00 (injected var(u_fov) = 0.04580)
  [CHECK] PASS  tau^2 is within the sampling band implied by the FOV count
            REML 0.04084 vs injected 0.04580 (ratio 0.89)
  [CHECK] PASS  Acquisition drift recovered with the injected sign
            realised slope -0.02331/FOV, recovered -0.02190
  [CHECK] PASS  HVG selection is enriched for the planted markers
            26/40 markers in the top 0.25 (50 genes); p = 1.166e-09
  [CHECK] PASS  Annotation recovers the simulated cell types
            ARI = 1.000 over 2000 cells
  ...
  13 passed, 0 failed, 2 informational, 0 skipped
```

Exact numbers shift slightly with BLAS threading and package versions; the
checks are written to tolerate that.

### Why the ground truth matters

The simulation is not arbitrary. It plants a per-FOV library-size offset of
known variance, a systematic decline in that offset across the acquisition
order, 40 marker genes carrying the between-cell-type variance, and five cell
types, then records the realised values in `ground_truth.json`. A demo run is
therefore checked for **correctness**, not merely for completing.
[`examples/demo-dataset/README.md`](examples/demo-dataset/README.md) documents
each check and the statistical reasoning behind its pass criterion.

The demo also doubles as the executable specification of the input schema in
the next section.

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
implementation of everything below. When this document and that script
disagree, the script is correct.

### 1. `filtered.zarr` schema

A [SpatialData](https://spatialdata.scverse.org) zarr archive. Only the table
is read; images, shapes and points are ignored (shape centroids are consulted
only as a last-resort fallback for coordinates).

**Structure**

| Element | Required | Notes |
|---|---|---|
| `tables['filtered']` | **yes** | the AnnData that the whole pipeline operates on |
| `tables['table']` | no | if present, re-indexed alongside `filtered` for MERSCOPE `EntityID` data |

**`X`: the expression matrix**

* **Raw integer counts**, cells × genes. Dense or `scipy.sparse`; both are handled.
* **Must not be normalized or log-transformed.** Every method starts from these
  counts, and TMM, scran and DESeq2 estimate size factors that are only
  meaningful on raw counts.
* Negative-control probes should already be **removed** from `X`; their
  per-cell total is carried in `obs` instead (see `total_counts_Negative`).

**`var`: genes**

* `var_names` holds gene symbols and must be unique.
* Symbols must overlap the scRNA-seq reference's row names, or annotation
  produces no usable genes and the run stops.
* The R stages normalize `-`, `_`, `:`, `/` and spaces to `.`
  (`clean_gene_names()` in `bin/annotation/run_insitutype.R`, and the same
  substitutions in `bin/plotters/hvg_upset_plot.py`). Symbols free of those
  characters are safest.

**`obs`: the index**

Cell IDs, unique and stable. They are the join key between the expression
matrix, the enriched metadata, the tissue annotation and every annotation
output, so they must survive a CSV round-trip unchanged. If an `EntityID`
column is present (the MERSCOPE convention) it is promoted to the index
automatically.

**`obs`: required columns**

All coordinates are in **micrometers**; there is no unit conversion anywhere in
the pipeline.

| Column | Type | Used by |
|---|---|---|
| `x_local_um`, `y_local_um` | float | pseudo-FOV rasterization; SpaNorm coordinates on Xenium/MERSCOPE |
| `x_global_um`, `y_global_um` | float | SpaNorm coordinates on CosMx; GeoJSON point-in-polygon; spatial plots |
| `area_um2` | float | `areaNorm` size factors. This method has no toggle, so the column is always required |

**`obs`: conditionally required**

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
| `.RData` | an object named exactly **`profile_matrix`**, genes × cell types, mean expression on the linear scale |
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

## Instructions for use

### Running BISTRO on your own data

**1. Prepare the inputs.** You need three files, described in full under
[Input data](#input-data):

- a quality-controlled `filtered.zarr`. BISTRO does no filtering of its own,
  so cell and gene QC must already be done
- an scRNA-seq reference (`.RData` holding `profile_matrix`, or a Seurat
  `.rds`), unless you set `skip_annotation = true`
- a tissue annotation (`.csv` or QuPath `.geojson`)

The single most common failure is a zarr that does not carry the required
`obs` columns. Check yours before launching a long run:

```python
import spatialdata as sd
obs = sd.SpatialData.read("my_filtered.zarr").tables["filtered"].obs
required = ["x_local_um", "y_local_um", "x_global_um", "y_global_um", "area_um2"]
print("missing:", [c for c in required if c not in obs.columns])
# plus 'total_counts_Negative' if you will run annotation,
# and 'fov' for platforms with native FOVs
```

[`examples/demo-dataset/make_demo_dataset.py`](examples/demo-dataset/make_demo_dataset.py)
builds a schema-correct archive from scratch and is the reference to copy
from.

**2. Copy the closest sample config and edit the paths.** Start from the entry
in [`sample_configs/`](sample_configs/) matching your platform, since
`technology`, `pixelSize`, `separate_fovs` and `add_global` are already set
appropriately there:

```bash
cp sample_configs/xenium-cancerBreast-5k.config my_dataset.config
```

Then edit, at minimum: `zarrFile`, `scReferenceFile`, `tissueAnnotation`,
`datasetName`, `output_folder_path`, `workDir` and `color_list`. Every
parameter is documented in [Configuration](#configuration).

**3. Size the resources.** The `process` block memory in the sample configs is
set for datasets of up to ~10⁶ cells × 18,000 genes. For a smaller panel,
reduce it; see the table under [Hardware](#hardware).

**4. Run.**

```bash
nextflow run run_BISTRO.nf -c my_dataset.config
```

On a cluster, add an executor and per-label environment setup. The pattern is
in [`examples/demo-dataset/demo.cluster.config`](examples/demo-dataset/demo.cluster.config),
which can be layered on with a second `-c` so the dataset config stays free of
site-specific settings:

```bash
nextflow run run_BISTRO.nf -c my_dataset.config -c my_cluster.config
```

**5. Read the report.** `output_folder_path/BISTRO_report.html` is
self-contained and collects every figure and table; see
[Outputs](#outputs) for the directory layout.

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

### Common problems

| Symptom | Cause |
|---|---|
| `No such variable: referenceAnnotationPath` | fixed; update to the current revision |
| `Unknown technology: ...` | `technology` must be exactly `CosMx`, `Xenium` or `MERSCOPE` |
| R package "not found" that is definitely installed | a missing system library, not a missing package. Run `Rscript envs/install_R_packages.R --check` |
| `Importing the numpy C-extensions failed` | `PYTHONPATH` from environment modules is shadowing the conda env; clear it (`env -u PYTHONPATH ...`) |
| Annotation step fails on missing `total_counts_Negative` | add the column, or set `skip_annotation = true` |
| `module load CMake` not found | only the pathway step does this; set `skip_pathway = true` off-cluster |

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

**`restore_published = true`** (default): each process checks
   `params.output_folder_path/<step>/` *before* running and short-circuits via
   symlinks if the expected output already exists. Set `restore_published =
   false` in the config to force a clean recomputation.

To rerun a single step from scratch, delete its sub-directory under the output
folder.

---

## Reproducing the manuscript results

Every quantitative result in the manuscript comes from running this pipeline
once per dataset, with the config committed here. There are no manual steps
between the pipeline output and the reported numbers.

> **Data availability.** The processed `filtered.zarr` archives, scRNA-seq
> references and tissue annotations are deposited at
> **[GEO accession, to be added on submission]**. The paths in
> `sample_configs/` point at the internal locations used during development;
> repoint `zarrFile`, `scReferenceFile`, `tissueAnnotation`,
> `output_folder_path` and `workDir` at your copy of the deposited data.

### The 20 datasets

| Config | `datasetName` | Platform | Panel | Tissue |
|---|---|---|---|---|
| `cosmx-crc-18k-11.config` | CosMx-ColonCancer-11-18k | CosMx | 18k | Colon |
| `cosmx-crc-18k-12.config` | CosMx-ColonCancer-12-18k | CosMx | 18k | Colon |
| `cosmx-crc-18k-21.config` | CosMx-ColonCancer-21-18k | CosMx | 18k | Colon |
| `cosmx-crc-18k-22.config` | CosMx-ColonCancer-22-18k | CosMx | 18k | Colon |
| `cosmx-crc-18k-23.config` | CosMx-ColonCancer-23-18k | CosMx | 18k | Colon |
| `cosmx-crc-18k-24.config` | CosMx-ColonCancer-24-18k | CosMx | 18k | Colon |
| `cosmx-cancerLiver-1k.config` | CosMx-LiverCancer-1k | CosMx | 1k | Liver (cancer) |
| `cosmx-normalLiver-1k.config` | CosMx-NormalLiver-1k | CosMx | 1k | Liver (normal) |
| `cosmx-normalPancreas-18k.config` | CosMx-Pancreas-18k | CosMx | 18k | Pancreas |
| `cosmx-normalPrefrontalCortex.config` | CosMx-PrefrontallCortex-6k | CosMx | 6k | Prefrontal cortex |
| `TMA_1404.config` | TMA1404 | CosMx | n/a | CRC liver-metastasis TMA |
| `merscope-brain-1k.config` | MERSCOPE-Brain-1k | MERSCOPE | 1k | Brain |
| `merscope-coloncancer-500-p1.config` | MERSCOPE-ColonCancer-P1-500 | MERSCOPE | 500 | Colon |
| `merscope-livercancer-500-p1.config` | MERSCOPE-LiverCancer-P1-500 | MERSCOPE | 500 | Liver |
| `merscope-lungcancer-500-p1.config` | MERSCOPE-LungCancer-P1-500 | MERSCOPE | 500 | Lung |
| `xenium-cancerBreast-5k.config` | Xenium-CancerBreast-5k | Xenium | 5k | Breast |
| `xenium-cancerBreastS1R1.config` | Xenium-CancerBreast-S1R1-300 | Xenium | 300 | Breast |
| `xenium-cancerBreastS1R2.config` | Xenium-CancerBreast-S1R2-300 | Xenium | 300 | Breast |
| `xenium-cancerBreastS2.config` | Xenium-CancerBreast-S2-300 | Xenium | 300 | Breast |
| `xenium-normalColon-300.config` | Xenium-HealthyColon-300 | Xenium | 300 | Colon |

`hvgThreshold = '0.75'` throughout. CosMx datasets use
`separate_fovs = '1'` and `add_global = '1'`; MERSCOPE and Xenium use `'0'`
and `'0'`. Three configs deviate deliberately:

- **`xenium-cancerBreast-5k.config`**: `separate_fovs = '1'` with
  `add_global = '0'`, the only Xenium dataset run per-FOV, because of the
  5,000-plex panel size.
- **`TMA_1404.config`**: the only one with `skip_hvg_bench = true` and with
  `vc_column = 'patient'`, since a TMA carries multiple patients per slide.
- **`cosmx-cancerLiver-1k.config`**: the only one setting
  `restore_published = true` explicitly.

### scRNA-seq references

| Reference | Used by |
|---|---|
| `ColonCRC.RData` | the six CosMx colon datasets, MERSCOPE-ColonCancer, Xenium-HealthyColon |
| `BreastCancer_Wu.RData` | all four Xenium breast datasets |
| `HCC_reference.RData` | CosMx-LiverCancer-1k, MERSCOPE-LiverCancer |
| `Brain_AllenBrainAtlas.RData` | CosMx-PrefrontalCortex, MERSCOPE-Brain |
| `Liver_HCA.RData` | CosMx-NormalLiver-1k |
| `Pancreas_HCA.RData` | CosMx-Pancreas-18k |
| `LungCancer.RData` | MERSCOPE-LungCancer |
| `LiuX_CRC_LiverMet_Reference_Profile.RData` | TMA1404 |

### Running them

One invocation per dataset:

```bash
for cfg in sample_configs/*.config; do
    nextflow run run_BISTRO.nf -c "$cfg"
done
```

In practice these were submitted individually, because a single 18k-plex dataset
needs up to 768 GB for the scran and SpaNorm steps and produces 100–400 GB of
intermediate matrices, so running all 20 concurrently is not realistic. Budget
several hours to a day per large dataset.

### Where each result comes from

| Manuscript section | Output |
|---|---|
| §2.2 FOV batch effect (τ², LRT, drift) | `evaluation/batch_effect/<datasetName>_batch_effect_summary.csv` |
| §2.3 Transformation / mean–variance | `evaluation/transformation/<datasetName>_transformation_summary.csv` |
| §2.4 HVG selection benchmark | `evaluation/hvg_benchmark/` |
| §2.4 Annotation agreement | `evaluation/annotation_agreement/<datasetName>_pairwise_ari_0.75.csv` |
| All figures | `BISTRO_report.html`, plus the per-step `plots/` directory |

Cross-dataset summary figures aggregate the per-dataset CSVs above over all 20
runs.

---

## License

Released under the MIT License. See [`LICENSE`](LICENSE).

BISTRO orchestrates third-party normalization methods that carry their own
licenses (among them edgeR, DESeq2, scran, SpaNorm, Seurat/SCTransform and
InSituType). Each is invoked as a separate process rather than linked, but if
you redistribute BISTRO together with those dependencies, check their terms.

---
