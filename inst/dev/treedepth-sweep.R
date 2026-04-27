#!/usr/bin/env Rscript
# Seed sweep harness for the default `estimate_infections_options` def chunk.
#
# Usage: Rscript inst/dev/treedepth-sweep.R LABEL ALPHA_SPEC SEED1 [SEED2 ...]
#
# ALPHA_SPEC is one of:
#   default                — keep gp_opts() defaults (Normal(0, 0.01))
#   normal:<sd>            — alpha = Normal(0, sd)
#   lognormal:<mean>:<sd>  — alpha = LogNormal(mean = mean, sd = sd) on natural scale
#
# Output: inst/dev/treedepth-sweep-results.csv (appended).

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(rstan)
  library(data.table)
})

options(mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript treedepth-sweep.R LABEL ALPHA_SPEC SEED1 [SEED2 ...]")
}
label <- args[[1]]
alpha_spec <- args[[2]]
seeds <- as.integer(args[-c(1, 2)])

parse_alpha <- function(spec) {
  if (spec == "default") return(NULL)
  parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
  family <- parts[[1]]
  if (family == "normal") {
    Normal(mean = 0, sd = as.numeric(parts[[2]]))
  } else if (family == "lognormal") {
    LogNormal(mean = as.numeric(parts[[2]]), sd = as.numeric(parts[[3]]))
  } else {
    stop("unknown alpha spec: ", spec)
  }
}
alpha_dist <- parse_alpha(alpha_spec)
gp_args <- if (is.null(alpha_dist)) gp_opts() else gp_opts(alpha = alpha_dist)
cat(sprintf("[%s] alpha_spec=%s\n", label, alpha_spec))

# vignette `def` chunk setup — exactly as in
# vignettes/estimate_infections_options.Rmd.orig
incubation_period <- LogNormal(
  meanlog = Normal(1.6, 0.05),
  sdlog = Normal(0.5, 0.05),
  max = 14
)
reporting_delay <- LogNormal(meanlog = 0.5, sdlog = 0.5, max = 10)
delay <- incubation_period + reporting_delay
rt_prior <- LogNormal(mean = 2, sd = 1)

out_path <- "inst/dev/treedepth-sweep-results.csv"
new_file <- !file.exists(out_path)
log_path <- "inst/dev/treedepth-sweep.log"
log_msg <- function(...) {
  cat(format(Sys.time()), ...,  "\n", file = log_path, append = TRUE)
}
log_msg(sprintf("=== Run started: label=%s alpha_spec=%s seeds=[%s]",
                label, alpha_spec, paste(seeds, collapse = ",")))

run_one <- function(seed) {
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
        warmup = 200,
        samples = 1200,
        seed = seed,
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
      min_bulk_ess = NA_real_,
      max_chain_sample_s = NA_real_
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
  rhat <- summ[, "Rhat"]
  bulk <- summ[, "n_eff"]
  data.table(
    label = label, seed = seed, elapsed_s = elapsed, ok = TRUE, msg = "",
    n_chains = n_chains, n_stuck = n_stuck,
    max_td = max(td), td_hits = sum(td >= 12),
    divergences = div,
    max_rhat = max(rhat, na.rm = TRUE),
    min_bulk_ess = min(bulk, na.rm = TRUE),
    max_chain_sample_s = max(el_t[, "sample"])
  )
}

results <- rbindlist(lapply(seeds, function(s) {
  cat(sprintf("[%s] seed=%d ...\n", label, s))
  log_msg(sprintf("[%s] seed=%d START", label, s))
  r <- run_one(s)
  summary_line <- sprintf(
    "stuck=%d/4 td_hits=%s div=%s max_rhat=%s max_chain_sample=%ss elapsed=%.1fs",
    r$n_stuck,
    if (is.na(r$td_hits)) "NA" else as.character(r$td_hits),
    if (is.na(r$divergences)) "NA" else as.character(r$divergences),
    if (is.na(r$max_rhat)) "NA" else sprintf("%.3g", r$max_rhat),
    if (is.na(r$max_chain_sample_s)) "NA" else sprintf("%.1f", r$max_chain_sample_s),
    r$elapsed_s
  )
  cat("  ", summary_line, "\n", sep = "")
  log_msg(sprintf("[%s] seed=%d DONE: %s", label, s, summary_line))
  r
}))

fwrite(results, out_path, append = !new_file)
cat(sprintf("Wrote %d rows to %s\n", nrow(results), out_path))
