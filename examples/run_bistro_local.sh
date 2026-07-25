#!/usr/bin/env bash
# ============================================================================
# BISTRO sample local launcher
# ----------------------------------------------------------------------------
# Runs the BISTRO Nextflow pipeline on a single machine 
#
# Usage:
#   1. Copy this script next to your data:
#        cp examples/run_bistro_local.sh my_run.sh
#   2. Edit the four paths in the USER CONFIGURATION block below.
#   3. Run it:
#        bash my_run.sh
#
# Notes:
#   * BISTRO is heavy: scran, SpaNorm and the HVG benchmark each peak at
#     several hundred GB of RAM on multi-thousand-cell datasets. For large
#     panels you almost certainly want an HPC / cloud node, not a laptop.
#   * Each dataset should use its OWN  NXF_CACHE  and  WORK_DIR  to avoid
#     Nextflow lock conflicts when running multiple BISTRO jobs in parallel.
# ============================================================================

set -euo pipefail

# ───── USER CONFIGURATION ───────────────────────────────────────────────────
# Absolute path to the cloned BISTRO repository (containing run_BISTRO.nf).
BISTRO_DIR="/path/to/BISTRO"

# Dataset config file. Pick one from sample_configs/ or write your own.
CONFIG_PATH="${BISTRO_DIR}/sample_configs/xenium-cancerBreast-5k.config"

# Per-run scratch directory. Nextflow writes all task work-dirs under here.
# Use a fresh directory for each dataset.
WORK_DIR="${HOME}/bistro_runs/xenium-cancerBreast-5k"

# Per-run Nextflow cache (.nextflow/ history, plugin cache, lock files).
NXF_CACHE="${HOME}/.nextflow/cache/xenium-cancerBreast-5k"
# ────────────────────────────────────────────────────────────────────────────

mkdir -p "${WORK_DIR}" "${NXF_CACHE}"
cd "${WORK_DIR}"

export NXF_HOME="${NXF_CACHE}"
export NXF_WORK="${WORK_DIR}/work"

echo "BISTRO launcher"
echo "  repo      : ${BISTRO_DIR}"
echo "  config    : ${CONFIG_PATH}"
echo "  work dir  : ${WORK_DIR}"
echo "  nxf cache : ${NXF_CACHE}"
echo

# Run Nextflow.
nextflow run "${BISTRO_DIR}/run_BISTRO.nf" \
    -c       "${CONFIG_PATH}" 

NF_EXIT=$?
echo ">>> Nextflow exited with code ${NF_EXIT}"
exit ${NF_EXIT}
