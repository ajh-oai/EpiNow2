#!/bin/bash
#SBATCH --job-name=epinow-td
#SBATCH --output=epinow-td_%j.out
#SBATCH --error=epinow-td_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail
cd ~/EpiNow2-td-investigation

# Compile package once before running so the rt.stan change is picked up
Rscript -e 'pkgbuild::compile_dll(quiet = FALSE)'
Rscript inst/dev/treedepth-experiments.R
