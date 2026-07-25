#!/usr/bin/env bash
# ============================================================================
# BISTRO demo -- SLURM submission (University of Bern / UBELIX)
# ----------------------------------------------------------------------------
# Runs the whole demo inside a SINGLE scheduler allocation: generate the
# simulated inputs, build the InSituType reference, run the pipeline, then
# verify the outputs against the injected ground truth.
#
#     sbatch examples/demo-dataset/run_demo_slurm.sbatch
#
# Watch it with:
#     tail -f <repo>/examples/demo-dataset/demo_slurm-<jobid>.out
#
# Off-cluster, use run_demo.sh instead -- it needs no scheduler and no modules.
# ============================================================================
#SBATCH --job-name=bistro-demo
#SBATCH --partition=epyc2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
# --nodes is mandatory whenever --mem is given on this cluster.

set -euo pipefail

# --- Locate the repository --------------------------------------------------
# SLURM copies the script, so BASH_SOURCE is unreliable; use the submission
# directory, falling back to this file's location for a direct bash run.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    REPO_DIR="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
    REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi

# Tolerate submission from either the repo root or the demo directory.
if [ ! -f "${REPO_DIR}/run_BISTRO.nf" ] && [ -f "${REPO_DIR}/../../run_BISTRO.nf" ]; then
    REPO_DIR="$(cd "${REPO_DIR}/../.." && pwd)"
fi

DEMO_DIR="${REPO_DIR}/examples/demo-dataset"

if [ ! -f "${REPO_DIR}/run_BISTRO.nf" ]; then
    echo "ERROR: run_BISTRO.nf not found under ${REPO_DIR}." >&2
    echo "Submit from the repository root: sbatch examples/demo-dataset/run_demo_slurm.sbatch" >&2
    exit 1
fi

echo "============================================================"
echo "BISTRO demo"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)"
echo "  repo     : ${REPO_DIR}"
echo "  demo dir : ${DEMO_DIR}"
echo "  started  : $(date)"
echo "============================================================"

# --- 1. Generate the simulated inputs (Python) ------------------------------
echo
echo "==> [1/4] Generating the simulated dataset"
(
    HPC_WORKSPACE=igmp_slide_workspace module load Workspace
    module load Anaconda3
    eval "$(conda shell.bash hook)"
    source activate py311
    python "${DEMO_DIR}/make_demo_dataset.py" --output_dir "${DEMO_DIR}"
)

# --- 2. Build the InSituType reference (R) ----------------------------------
echo
echo "==> [2/4] Building the InSituType reference"
(
    HPC_WORKSPACE=igmp_slide_workspace module load Workspace
    HPC_WORKSPACE=igmp_dp_workspace module load Workspace
    module load GCCcore/12.3.0
    module load R/4.4.2-gfbf-2024a
    export R_LIBS_USER="/storage/research/igmp_dp_workspace/carreno_jose/R/x86_64-pc-linux-gnu-library/4.4"
    Rscript "${DEMO_DIR}/make_demo_reference.R" \
        "${DEMO_DIR}/demo_reference_profiles.csv" \
        "${DEMO_DIR}/demo_reference.RData"
)

# --- 3. Run the pipeline ----------------------------------------------------
# Nextflow itself only needs Java. Each process loads its own environment via
# the beforeScript blocks in demo.cluster.config, so R and conda never collide.
echo
echo "==> [3/4] Running BISTRO"
module load Java/17.0.6
module load Nextflow/24.04.4

cd "${REPO_DIR}"
nextflow run "${REPO_DIR}/run_BISTRO.nf" \
    -c "${DEMO_DIR}/demo.config" \
    -c "${DEMO_DIR}/demo.cluster.config" \
    -ansi-log false

# --- 4. Verify against the injected ground truth ----------------------------
echo
echo "==> [4/4] Verifying outputs against the ground truth"
(
    HPC_WORKSPACE=igmp_slide_workspace module load Workspace
    module load Anaconda3
    eval "$(conda shell.bash hook)"
    source activate py311
    python "${DEMO_DIR}/check_demo.py" \
        --demo_dir "${DEMO_DIR}" \
        --output_dir "${DEMO_DIR}/demo_output"
)
VERIFY_EXIT=$?

echo
echo "============================================================"
echo "  finished : $(date)"
echo "  report   : ${DEMO_DIR}/demo_output/BISTRO_report.html"
echo "  verify   : exit ${VERIFY_EXIT} (0 = all checks passed)"
echo "============================================================"

exit ${VERIFY_EXIT}
