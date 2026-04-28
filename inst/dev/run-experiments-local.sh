#!/bin/bash
# Run the array-style sweep on a single fast machine. Loops over (config x seed)
# sequentially within each pass; recompiles between baseline and rescaled rt.stan.
#
# Usage:
#   cd path/to/EpiNow2  (the branch investigate/gp-alpha-identifiability)
#   bash inst/dev/run-experiments-local.sh
#
# Env overrides: SEEDS_CSV (default "1,2,3,...,15"), CONFIGS_BASELINE,
# CONFIGS_RESCALED.
#
# Outputs:
#   inst/dev/array-results/${config}/${config}_seed${seed}.csv   one per fit
#   inst/dev/array-results/run.log                              stdout/stderr

set -euo pipefail

SEEDS_CSV=${SEEDS_CSV:-"1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"}
CONFIGS_BASELINE=${CONFIGS_BASELINE:-"baseline_default baseline_ln0p3"}
CONFIGS_RESCALED=${CONFIGS_RESCALED:-"rescaled_default rescaled_ln0p3"}

mkdir -p inst/dev/array-results
LOG=inst/dev/array-results/run.log
date | tee -a "$LOG"

RT_STAN=inst/stan/functions/rt.stan
BACKUP=inst/stan/functions/rt.stan.rescaled
cp $RT_STAN $BACKUP

force_recompile() {
  rm -f src/stanExports_estimate_infections.{cc,h,o} src/EpiNow2.so
  Rscript -e 'rstantools::rstan_config(".")' >/dev/null 2>&1
  Rscript -e 'pkgbuild::compile_dll(quiet = TRUE)' 2>&1 | tail -1
}

run_pass() {
  local pass_label=$1; shift
  local configs="$@"
  echo "=== Recompile for pass=$pass_label ===" | tee -a "$LOG"
  force_recompile 2>&1 | tee -a "$LOG"
  for cfg in $configs; do
    mkdir -p "inst/dev/array-results/$cfg"
    for seed in $(echo "$SEEDS_CSV" | tr ',' ' '); do
      out="inst/dev/array-results/$cfg/${cfg}_seed$(printf '%05d' $seed).csv"
      if [ -s "$out" ]; then
        echo "skip $cfg seed=$seed (csv exists)" | tee -a "$LOG"
        continue
      fi
      echo "$(date) starting $cfg seed=$seed" | tee -a "$LOG"
      SLURM_ARRAY_TASK_ID=$seed Rscript inst/dev/treedepth-array-fit.R \
        "$cfg" "inst/dev/array-results/$cfg" 2>&1 | tee -a "$LOG"
    done
  done
}

# Pass 1: baseline rt.stan
sed -i 's|gp\[2:(gp_n + 1)\] = noise / sqrt(gp_n);|gp[2:(gp_n + 1)] = noise;|' $RT_STAN
sed -i '/alpha controls trajectory SD/,/approximately alpha\^2 instead/d' $RT_STAN
run_pass baseline $CONFIGS_BASELINE

# Pass 2: rescaled rt.stan
cp $BACKUP $RT_STAN
run_pass rescaled $CONFIGS_RESCALED

echo "=== Done at $(date) ===" | tee -a "$LOG"
echo "Aggregate with: cat inst/dev/array-results/*/*.csv | sort -u" | tee -a "$LOG"
