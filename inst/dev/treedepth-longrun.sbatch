#!/bin/bash
#SBATCH --job-name=epinow-td-long
#SBATCH --output=epinow-td-long_%j.out
#SBATCH --error=epinow-td-long_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Run a single fit with NO max_execution_time cap and a generous wallclock to
# determine actual fit duration on HPC. If fits genuinely complete in under
# 10 min, then the previous 600s cap was wrong; if not, the model truly is
# pathological at all alpha priors tested.

set -euo pipefail
cd ~/EpiNow2-td-investigation
module load R/4.4.0

cat > /tmp/longrun.R <<'EOF'
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(rstan)
})
options(mc.cores = 4)
incubation_period <- LogNormal(meanlog = Normal(1.6, 0.05), sdlog = Normal(0.5, 0.05), max = 14)
reporting_delay <- LogNormal(meanlog = 0.5, sdlog = 0.5, max = 10)
delay <- incubation_period + reporting_delay
rt_prior <- LogNormal(mean = 2, sd = 1)

# rescaled GP + LogNormal alpha bounded away from 0 — the most promising
# combination from the experiments
gp <- gp_opts(alpha = LogNormal(mean = 0.3, sd = 0.3))

cat("Starting fit at:", format(Sys.time()), "\n"); t0 <- Sys.time()
fit_obj <- estimate_infections(
  example_confirmed,
  generation_time = gt_opts(example_generation_time),
  delays = delay_opts(delay),
  rt = rt_opts(prior = rt_prior),
  gp = gp,
  stan = stan_opts(warmup = 250, samples = 2000, seed = 42L,
                   control = list(adapt_delta = 0.9, max_treedepth = 12))
)
cat("Total elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
fit <- fit_obj$fit
cat("Per-chain elapsed:\n"); print(rstan::get_elapsed_time(fit))
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
td <- sapply(sp, function(x) x[, "treedepth__"])
cat("td>=12 per chain:", apply(td, 2, function(x) sum(x >= 12)), "\n")
cat("max td per chain:", apply(td, 2, max), "\n")
cat("div per chain:", sapply(sp, function(x) sum(x[, "divergent__"])), "\n")
EOF

Rscript /tmp/longrun.R
