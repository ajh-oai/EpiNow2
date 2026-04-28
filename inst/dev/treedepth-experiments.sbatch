#!/bin/bash
#SBATCH --job-name=epinow-td-r
#SBATCH --output=epinow-td-r_%j.out
#SBATCH --error=epinow-td-r_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Re-runs ONLY the rescaled pass with a forced full recompile. The previous
# version of this script reused stanExports_estimate_infections.o because
# pkgbuild::compile_dll does not re-trigger rstan_config when the only change
# is in an included .stan file.

set -euo pipefail
cd ~/EpiNow2-td-investigation

module load R/4.4.0

# Force full recompile: regenerate stanExports_*.{cc,h} from current .stan
# code via rstan_config, and remove all stale .o so they get rebuilt.
rm -f src/stanExports_estimate_infections.{cc,h,o}
rm -f src/EpiNow2.so

echo "=== Reconfiguring rstan models from current .stan code ==="
Rscript -e 'rstantools::rstan_config(".")'
ls -la src/stanExports_estimate_infections.{cc,h} || true

echo "=== Compiling rescaled rt.stan from clean ==="
Rscript -e 'pkgbuild::compile_dll(quiet = FALSE)' 2>&1 | tail -5
ls -la src/stanExports_estimate_infections.o src/EpiNow2.so

echo "=== Running rescaled experiments ==="
Rscript inst/dev/treedepth-experiments.R rescaled

echo "=== Done ==="
