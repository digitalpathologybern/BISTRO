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
# PYTHONPATH is cleared for every Python call below. On systems with
# environment modules, R's spatial dependencies (GDAL in particular) put their
# own site-packages on PYTHONPATH, and a conda interpreter honours it: the
# result is that Python imports a numpy built for a different interpreter and
# dies with "Importing the numpy C-extensions failed". Clearing it costs
# nothing on a normal desktop, where PYTHONPATH is usually unset anyway.
echo "==> [1/4] Generating the simulated dataset"
env -u PYTHONPATH python "${DEMO_DIR}/make_demo_dataset.py" --output_dir "${DEMO_DIR}"
echo

# ---- 2. Convert the reference profiles to the .RData InSituType expects -----
echo "==> [2/4] Building the InSituType reference"
Rscript "${DEMO_DIR}/make_demo_reference.R" \
    "${DEMO_DIR}/demo_reference_profiles.csv" \
    "${DEMO_DIR}/demo_reference.RData"
echo

# ---- 3. Run the pipeline ---------------------------------------------------
# PYTHONPATH is cleared here too, not just for the calls above: Nextflow hands
# its environment to every task it spawns, so a polluted PYTHONPATH would
# otherwise resurface inside the pipeline processes rather than in this script.
echo "==> [3/4] Running BISTRO"
cd "${REPO_DIR}"
env -u PYTHONPATH nextflow run "${REPO_DIR}/run_BISTRO.nf" -c "${DEMO_DIR}/demo.config"
echo

# ---- 4. Verify against the injected ground truth ----------------------------
echo "==> [4/4] Verifying outputs against the ground truth"
env -u PYTHONPATH python "${DEMO_DIR}/check_demo.py" \
    --demo_dir "${DEMO_DIR}" \
    --output_dir "${DEMO_DIR}/demo_output"
VERIFY_EXIT=$?

echo
echo "Demo complete."
echo "  report        : ${DEMO_DIR}/demo_output/BISTRO_report.html"
echo "  ground truth  : ${DEMO_DIR}/ground_truth.json"
echo "  verification  : exit ${VERIFY_EXIT} (0 = all checks passed)"
echo
echo "The checks are explained in ${DEMO_DIR}/README.md"

exit ${VERIFY_EXIT}
