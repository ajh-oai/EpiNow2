#!/usr/bin/env Rscript
# Single-fit runner used by the SLURM array job.
# One fit per SLURM array task — uses ${SLURM_ARRAY_TASK_ID} as seed and
# argv[1] to pick the experiment config.
#
# Usage: Rscript treedepth-array-fit.R CONFIG OUTDIR
#   CONFIG: baseline_default | baseline_ln0p3 | rescaled_default | rescaled_ln0p3
#   OUTDIR: directory to write per-task csv row to (parent collects)

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(rstan)
  library(data.table)
})
options(mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
config <- args[[1]]
outdir <- args[[2]]
seed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

alpha_dist <- switch(
  config,
  baseline_default = NULL,
  baseline_ln0p3 = LogNormal(mean = 0.3, sd = 0.3),
  rescaled_default = NULL,
  rescaled_ln0p3 = LogNormal(mean = 0.3, sd = 0.3),
  stop("unknown config: ", config)
)
gp_args <- if (is.null(alpha_dist)) gp_opts() else gp_opts(alpha = alpha_dist)

incubation_period <- LogNormal(meanlog = Normal(1.6, 0.05), sdlog = Normal(0.5, 0.05), max = 14)
reporting_delay <- LogNormal(meanlog = 0.5, sdlog = 0.5, max = 10)
delay <- incubation_period + reporting_delay
rt_prior <- LogNormal(mean = 2, sd = 1)

cat(sprintf("[%s seed=%d] START at %s\n", config, seed, format(Sys.time())))
t0 <- Sys.time()
fit_obj <- tryCatch(
  estimate_infections(
    example_confirmed,
    generation_time = gt_opts(example_generation_time),
    delays = delay_opts(delay),
    rt = rt_opts(prior = rt_prior),
    gp = gp_args,
    stan = stan_opts(
      warmup = 250, samples = 2000, seed = seed,
      control = list(adapt_delta = 0.9, max_treedepth = 11)
    )
  ),
  error = function(e) e
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if (inherits(fit_obj, "error")) {
  row <- data.table(
    config = config, seed = seed, elapsed_s = elapsed, ok = FALSE,
    msg = conditionMessage(fit_obj),
    n_chains = NA_integer_, n_stuck = NA_integer_,
    td_hits = NA_integer_, max_td = NA_integer_,
    divergences = NA_integer_, max_rhat = NA_real_,
    min_bulk_ess = NA_real_,
    max_chain_warmup_s = NA_real_, max_chain_sample_s = NA_real_
  )
} else {
  fit <- fit_obj$fit
  n_chains <- fit@sim$chains
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  td <- unlist(lapply(sp, function(x) x[, "treedepth__"]))
  div <- sum(unlist(lapply(sp, function(x) x[, "divergent__"])))
  el_t <- rstan::get_elapsed_time(fit)
  summ <- rstan::summary(fit)$summary
  row <- data.table(
    config = config, seed = seed, elapsed_s = elapsed, ok = TRUE, msg = "",
    n_chains = n_chains, n_stuck = 4L - n_chains,
    td_hits = sum(td >= 11), max_td = max(td),
    divergences = div,
    max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
    min_bulk_ess = min(summ[, "n_eff"], na.rm = TRUE),
    max_chain_warmup_s = max(el_t[, "warmup"]),
    max_chain_sample_s = max(el_t[, "sample"])
  )
}

out_path <- file.path(outdir, sprintf("%s_seed%05d.csv", config, seed))
fwrite(row, out_path)
cat(sprintf("[%s seed=%d] DONE: ok=%s elapsed=%.0fs td_hits=%s div=%s max_chain_sample=%s\n",
            config, seed, row$ok,
            row$elapsed_s,
            ifelse(is.na(row$td_hits), "NA", as.character(row$td_hits)),
            ifelse(is.na(row$divergences), "NA", as.character(row$divergences)),
            ifelse(is.na(row$max_chain_sample_s), "NA", sprintf("%.1f", row$max_chain_sample_s))))
