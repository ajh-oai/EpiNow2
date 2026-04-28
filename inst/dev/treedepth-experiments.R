#!/usr/bin/env Rscript
# Run baseline + all experiments in one session.
# Each experiment runs N seeds, recording per-chain treedepth and divergences.
# Output: inst/dev/treedepth-experiments-results.csv
#
# Designed to be run on a faster machine (e.g. HPC) where Stan model fit times
# are short enough to do meaningful seed sweeps.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(rstan)
  library(data.table)
})

options(mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[[1]] else "baseline"

seeds <- c(42L, 9876L, 100L, 200L, 300L)

# Sweep over alpha priors. Same set is run for both modes ("baseline" and
# "rescaled") so cross-mode comparison at fixed alpha is meaningful.
experiments <- list(
  list(label = "alpha_default",  alpha_dist = NULL),
  list(label = "alpha_n0p05",    alpha_dist = Normal(mean = 0, sd = 0.05)),
  list(label = "alpha_n0p1",     alpha_dist = Normal(mean = 0, sd = 0.1)),
  list(label = "alpha_n0p3",     alpha_dist = Normal(mean = 0, sd = 0.3)),
  list(label = "alpha_ln0p3",    alpha_dist = LogNormal(mean = 0.3, sd = 0.3))
)

# vignette `def` chunk setup
incubation_period <- LogNormal(
  meanlog = Normal(1.6, 0.05), sdlog = Normal(0.5, 0.05), max = 14
)
reporting_delay <- LogNormal(meanlog = 0.5, sdlog = 0.5, max = 10)
delay <- incubation_period + reporting_delay
rt_prior <- LogNormal(mean = 2, sd = 1)

out_path <- sprintf("inst/dev/treedepth-experiments-%s.csv", mode)
log_path <- sprintf("inst/dev/treedepth-experiments-%s.log", mode)
log_msg <- function(...) {
  cat(format(Sys.time()), ..., "\n", file = log_path, append = TRUE)
}

run_one <- function(label, alpha_dist, seed) {
  gp_args <- if (is.null(alpha_dist)) gp_opts() else gp_opts(alpha = alpha_dist)
  set.seed(seed)
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
        max_execution_time = 600,
        control = list(adapt_delta = 0.9, max_treedepth = 12)
      )
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit_obj, "error")) {
    return(data.table(
      label = label, seed = seed, elapsed_s = elapsed,
      ok = FALSE, msg = conditionMessage(fit_obj),
      n_chains = NA_integer_, n_stuck = 4L,
      max_td = NA_integer_, td_hits = NA_integer_,
      divergences = NA_integer_, max_rhat = NA_real_,
      min_bulk_ess = NA_real_, max_chain_sample_s = NA_real_
    ))
  }
  fit <- fit_obj$fit
  n_chains <- fit@sim$chains
  n_stuck <- 4L - n_chains
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  td <- unlist(lapply(sp, function(x) x[, "treedepth__"]))
  div <- sum(unlist(lapply(sp, function(x) x[, "divergent__"])))
  el_t <- rstan::get_elapsed_time(fit)
  summ <- rstan::summary(fit)$summary
  data.table(
    label = label, seed = seed, elapsed_s = elapsed, ok = TRUE, msg = "",
    n_chains = n_chains, n_stuck = n_stuck,
    max_td = max(td), td_hits = sum(td >= 12),
    divergences = div,
    max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
    min_bulk_ess = min(summ[, "n_eff"], na.rm = TRUE),
    max_chain_sample_s = max(el_t[, "sample"])
  )
}

results <- list()
for (exp in experiments) {
  for (s in seeds) {
    full_label <- sprintf("%s__%s", mode, exp$label)
    log_msg(sprintf("[%s] seed=%d START", full_label, s))
    r <- run_one(full_label, exp$alpha_dist, s)
    log_msg(sprintf(
      "[%s] seed=%d DONE: stuck=%d/4 td_hits=%s div=%s max_rhat=%s elapsed=%.0fs",
      full_label, s, r$n_stuck,
      if (is.na(r$td_hits)) "NA" else as.character(r$td_hits),
      if (is.na(r$divergences)) "NA" else as.character(r$divergences),
      if (is.na(r$max_rhat)) "NA" else sprintf("%.3g", r$max_rhat),
      r$elapsed_s
    ))
    results[[length(results) + 1]] <- r
    fwrite(rbindlist(results), out_path)  # incremental write
  }
}

cat("Wrote", nrow(rbindlist(results)), "rows to", out_path, "\n")
