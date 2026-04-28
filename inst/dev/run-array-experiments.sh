#!/bin/bash
# Submit two array passes (baseline rt.stan vs rescaled rt.stan), each with
# 4 alpha-prior configs × 15 seeds = 60 fits per pass, all running as SLURM
# array tasks in parallel. Total = 120 array tasks. Each task gets 4 cpus +
# 2h budget; SLURM may schedule them across many nodes simultaneously.
#
# Run once for baseline, switch rt.stan, run again for rescaled.

set -euo pipefail
cd ~/EpiNow2-td-investigation
module load R/4.4.0

RT_STAN=inst/stan/functions/rt.stan
BACKUP=inst/stan/functions/rt.stan.rescaled
cp $RT_STAN $BACKUP

force_recompile() {
  rm -f src/stanExports_estimate_infections.{cc,h,o} src/EpiNow2.so
  Rscript -e 'rstantools::rstan_config(".")' >/dev/null 2>&1
  Rscript -e 'pkgbuild::compile_dll(quiet = TRUE)' 2>&1 | tail -1
}

submit_pass() {
  local pass_label=$1
  for cfg in ${pass_label}_default ${pass_label}_ln0p3; do
    sbatch --job-name="td-${cfg}" --export=CONFIG=$cfg \
      inst/dev/treedepth-array.sbatch | tee -a inst/dev/array-jobs.log
  done
}

# ===== PASS 1: baseline (no rescale) =====
sed -i 's|gp\[2:(gp_n + 1)\] = noise / sqrt(gp_n);|gp[2:(gp_n + 1)] = noise;|' $RT_STAN
sed -i '/alpha controls trajectory SD/,/approximately alpha\^2 instead/d' $RT_STAN
force_recompile
echo "=== Submitting baseline arrays ==="
submit_pass baseline

# Wait for the baseline arrays to complete before recompiling, otherwise
# tasks scheduled later would link the rescaled .so.
echo "=== Waiting for baseline tasks to finish ==="
while squeue -u "$USER" -h -o "%j" | grep -q "td-baseline"; do
  sleep 60
done
echo "=== Baseline pass complete ==="

# ===== PASS 2: rescaled =====
cp $BACKUP $RT_STAN
force_recompile
echo "=== Submitting rescaled arrays ==="
submit_pass rescaled

echo "=== All arrays submitted ==="
