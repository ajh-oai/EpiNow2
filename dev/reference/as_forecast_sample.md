# Convert EpiNow2 model output to a `forecast_sample` object

**\[experimental\]** Convert outputs of EpiNow2 fitting and forecasting
functions to `forecast_sample` objects via
[`scoringutils::as_forecast_sample()`](https://epiforecasts.io/scoringutils/reference/as_forecast_sample.html)
for evaluating predictive performance. Methods are provided for objects
returned by
[`epinow()`](https://epiforecasts.io/EpiNow2/dev/reference/epinow.md),
[`estimate_infections()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_infections.md),
[`forecast_secondary()`](https://epiforecasts.io/EpiNow2/dev/reference/forecast_secondary.md),
and
[`estimate_truncation()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_truncation.md).

These methods extract sample-level posterior predictions via
[`get_predictions()`](https://epiforecasts.io/EpiNow2/dev/reference/get_predictions.md)
with `format = "sample"`, merge them with the supplied observations on
`date`, and pass the result to
[`scoringutils::as_forecast_sample()`](https://epiforecasts.io/scoringutils/reference/as_forecast_sample.html).

[scoringutils](https://epiforecasts.io/scoringutils/reference/scoringutils-package.html)
is an optional dependency; calling these methods without it installed
gives an informative error.

## Usage

``` r
# S3 method for class 'estimate_infections'
as_forecast_sample(data, observations, horizon = 0, ...)

# S3 method for class 'epinow'
as_forecast_sample(data, observations, horizon = 0, ...)

# S3 method for class 'forecast_secondary'
as_forecast_sample(data, observations, horizon = 0, ...)

# S3 method for class 'estimate_truncation'
as_forecast_sample(data, observations, horizon = -Inf, ...)
```

## Arguments

- data:

  Output of
  [`epinow()`](https://epiforecasts.io/EpiNow2/dev/reference/epinow.md),
  [`estimate_infections()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_infections.md),
  [`forecast_secondary()`](https://epiforecasts.io/EpiNow2/dev/reference/forecast_secondary.md),
  or
  [`estimate_truncation()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_truncation.md).

- observations:

  A `<data.frame>` of observed values to score against. Must contain a
  `date` column. For
  [`epinow()`](https://epiforecasts.io/EpiNow2/dev/reference/epinow.md)
  and
  [`estimate_infections()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_infections.md)
  objects must also contain a `confirm` column; for
  [`forecast_secondary()`](https://epiforecasts.io/EpiNow2/dev/reference/forecast_secondary.md)
  objects a `secondary` column; for
  [`estimate_truncation()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_truncation.md)
  objects a `confirm` column representing the latest, least-truncated
  observations.

- horizon:

  Numeric scalar lower bound on the `horizon` column of
  [`get_predictions()`](https://epiforecasts.io/EpiNow2/dev/reference/get_predictions.md)
  output. Predictions with a `horizon` value at or above this bound are
  retained. Defaults to `0` for
  [`epinow()`](https://epiforecasts.io/EpiNow2/dev/reference/epinow.md),
  [`estimate_infections()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_infections.md)
  and
  [`forecast_secondary()`](https://epiforecasts.io/EpiNow2/dev/reference/forecast_secondary.md)
  (i.e. forecast period only) and to `-Inf` for
  [`estimate_truncation()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_truncation.md)
  (keep all reconstructed horizons). Pass `horizon = -Inf` to disable
  filtering.

- ...:

  Additional arguments passed to
  [`scoringutils::as_forecast_sample()`](https://epiforecasts.io/scoringutils/reference/as_forecast_sample.html).
  `forecast_unit` is set automatically from the object class
  (`forecast_date`, `date`, `horizon`, plus `dataset` for
  [`estimate_truncation()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_truncation.md))
  and cannot be overridden.

## Value

A `forecast_sample` object as returned by
[`scoringutils::as_forecast_sample()`](https://epiforecasts.io/scoringutils/reference/as_forecast_sample.html).
Rows for which `observations` does not provide a value on the
corresponding `date` are dropped.

## See also

[`get_predictions()`](https://epiforecasts.io/EpiNow2/dev/reference/get_predictions.md)
for the underlying sample extraction.

## Examples

``` r
# \donttest{
library(scoringutils)

# samples and calculation time have been reduced for this example
# for real analyses, use at least samples = 2000
fit <- estimate_infections(example_confirmed[1:40],
  generation_time = gt_opts(example_generation_time),
  delays = delay_opts(example_incubation_period + example_reporting_delay),
  rt = rt_opts(prior = LogNormal(mean = 2, sd = 0.2)),
  stan = stan_opts(samples = 100, warmup = 200)
)
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

forecast_obj <- as_forecast_sample(fit, observations = example_confirmed)
score(forecast_obj)
#> Warning: Predictions appear to be integer-valued.
#> ! The log score uses kernel density estimation, which may not be appropriate
#>   for integer-valued forecasts.
#> ℹ See the scoringRules package for alternatives for discrete probability
#>   distributions.
#>    forecast_date       date horizon  bias      dss     crps overprediction
#>           <Date>     <Date>   <num> <num>    <num>    <num>          <num>
#> 1:    2020-04-01 2020-04-01       0 -0.02 13.78257 232.7184           0.00
#> 2:    2020-04-01 2020-04-02       1 -0.28 14.15587 337.7138           0.00
#> 3:    2020-04-01 2020-04-03       2  0.32 15.03383 468.8191          80.14
#> 4:    2020-04-01 2020-04-04       3 -0.30 14.35096 365.9270           0.00
#> 5:    2020-04-01 2020-04-05       4 -0.30 14.49760 351.7080           0.00
#> 6:    2020-04-01 2020-04-06       5 -0.08 14.59981 316.7176           0.00
#> 7:    2020-04-01 2020-04-07       6 -0.08 14.06452 235.9707           0.00
#> 8:    2020-04-01 2020-04-08       7 -0.26 13.79265 289.6675           0.00
#>    underprediction dispersion log_score       mad ae_median     se_mean
#>              <num>      <num>     <num>     <num>     <num>       <num>
#> 1:            0.28   232.4384  7.894147  991.1181      30.0    171.0864
#> 2:           58.68   279.0338  8.159316 1187.5626     465.0 149567.8276
#> 3:            0.00   388.6791  8.408873 1729.4529     483.5 560776.3225
#> 4:           87.00   278.9270  8.219774 1201.6473     478.5 124679.6100
#> 5:           64.18   287.5280  8.169082 1158.6519     417.0  97956.4804
#> 6:            5.70   311.0176  8.197606 1257.2448     105.0  12683.2644
#> 7:            5.30   230.6707  7.856085  985.9290      81.0   8029.9521
#> 8:           71.32   218.3475  8.070536  814.6887     419.5  18887.0049
# }
```
