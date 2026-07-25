#!/usr/bin/env bash
# ============================================================================
# BISTRO demo runner
# ----------------------------------------------------------------------------
# Generates the simulated demo dataset and runs the full BISTRO pipeline on it,
# locally, in minutes.
#
#   bash examples/demo-dataset/run_demo.sh
#
# Prerequisites: the conda environment from envs/environment.yml active, the R
# packages from envs/R_packages.csv installed, and nextflow on PATH. See the
# "Demo" section of the top-level README.
#
# Results land in examples/demo-dataset/demo_output/, with the HTML report at
# examples/demo-dataset/demo_output/BISTRO_report.html.
# ============================================================================

set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${DEMO_DIR}/../.." && pwd)"

echo "BISTRO demo"
echo "  repo      : ${REPO_DIR}"
echo "  demo dir  : ${DEMO_DIR}"
echo

# ---- Check prerequisites ---------------------------------------------------
missing=0
for cmd in python Rscript nextflow; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "ERROR: '${cmd}' not found on PATH." >&2
        missing=1
    fi
done
if [ "${missing}" -ne 0 ]; then
    echo >&2
    echo "Activate the BISTRO environment before running the demo:" >&2
    echo "    conda env create -n bistro -f ${REPO_DIR}/envs/environment.yml" >&2
    echo "    conda activate bistro" >&2
    exit 1
fi

# ---- 1. Generate the simulated inputs --------------------------------------
echo "==> [1/3] Generating the simulated dataset"
python "${DEMO_DIR}/make_demo_dataset.py" --output_dir "${DEMO_DIR}"
echo

# ---- 2. Convert the reference profiles to the .RData InSituType expects -----
echo "==> [2/3] Building the InSituType reference"
Rscript "${DEMO_DIR}/make_demo_reference.R" \
    "${DEMO_DIR}/demo_reference_profiles.csv" \
    "${DEMO_DIR}/demo_reference.RData"
echo

# ---- 3. Run the pipeline ---------------------------------------------------
echo "==> [3/3] Running BISTRO"
cd "${REPO_DIR}"
nextflow run "${REPO_DIR}/run_BISTRO.nf" -c "${DEMO_DIR}/demo.config"

echo
echo "Demo complete."
echo "  report        : ${DEMO_DIR}/demo_output/BISTRO_report.html"
echo "  ground truth  : ${DEMO_DIR}/ground_truth.json"
echo
echo "Compare the run against the injected ground truth as described in"
echo "  ${DEMO_DIR}/README.md"
