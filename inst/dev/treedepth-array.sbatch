#!/bin/bash
#SBATCH --job-name=epinow-arr
#SBATCH --output=arr_%A_%a.out
#SBATCH --error=arr_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-15

# Array job: one fit per task, seed = SLURM_ARRAY_TASK_ID. Caller sets
# CONFIG environment variable before sbatch (e.g. CONFIG=baseline_default).
# Per-task wallclock 2h covers even a stuck-chain fit.

set -euo pipefail
cd ~/EpiNow2-td-investigation
module load R/4.4.0

CONFIG=${CONFIG:?CONFIG must be set in environment}
OUTDIR=inst/dev/array-results/${CONFIG}
mkdir -p "$OUTDIR"

Rscript inst/dev/treedepth-array-fit.R "$CONFIG" "$OUTDIR"
