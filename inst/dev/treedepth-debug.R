#!/usr/bin/env Rscript
# Debug: run a single fit with verbose output.
#
# Usage: Rscript inst/dev/treedepth-debug.R LABEL ALPHA_SPEC SEED
#   ALPHA_SPEC: default | normal:<sd> | lognormal:<mean>:<sd>

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(rstan)
})

options(mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1) args[[1]] else "debug"
alpha_spec <- if (length(args) >= 2) args[[2]] else "default"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 42L

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

incubation_period <- LogNormal(
  meanlog = Normal(1.6, 0.05),
  sdlog = Normal(0.5, 0.05),
  max = 14
)
reporting_delay <- LogNormal(meanlog = 0.5, sdlog = 0.5, max = 10)
delay <- incubation_period + reporting_delay
rt_prior <- LogNormal(mean = 2, sd = 1)

cat(sprintf("[%s] seed=%d alpha_spec=%s\n", label, seed, alpha_spec))
t0 <- Sys.time()
fit_obj <- estimate_infections(
  example_confirmed,
  generation_time = gt_opts(example_generation_time),
  delays = delay_opts(delay),
  rt = rt_opts(prior = rt_prior),
  gp = gp_args,
  stan = stan_opts(
    warmup = 250,
    samples = 2000,
    seed = seed
  )
)
cat("Total elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
fit <- fit_obj$fit
cat("Per-chain elapsed (warmup, sample):\n")
print(rstan::get_elapsed_time(fit))

sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
td <- sapply(sp, function(x) x[, "treedepth__"])
cat("Treedepth >= 12 hits per chain:", apply(td, 2, function(x) sum(x >= 12)), "\n")
cat("Max treedepth per chain:", apply(td, 2, max), "\n")
div <- sapply(sp, function(x) sum(x[, "divergent__"]))
cat("Divergences per chain:", div, "\n")

# Look at posterior of alpha and rho across chains
alpha_post <- rstan::extract(fit, pars = "params", permuted = FALSE)
cat("dim(alpha_post):\n"); print(dim(alpha_post))
# try to look up alpha index by name
pnames <- dimnames(alpha_post)$parameters
cat("Parameter names sample:\n"); print(head(pnames, 10))
