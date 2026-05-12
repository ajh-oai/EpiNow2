# Probability distributions

Probability distributions

Generates a nonparametric distribution.

## Usage

``` r
LogNormal(meanlog, sdlog, mean, sd, ...)

Gamma(shape, rate, scale, mean, sd, ...)

Normal(mean, sd, ...)

Exp(rate, mean, ...)

Weibull(shape, scale, mean, sd, ...)

NonParametric(pmf, ...)

Fixed(value, ...)
```

## Arguments

- meanlog, sdlog:

  mean and standard deviation of the distribution on the log scale with
  default values of `0` and `1` respectively.

- mean, sd:

  mean and standard deviation of the distribution

- ...:

  arguments to define the limits of the distribution that will be passed
  to
  [`bound_dist()`](https://epiforecasts.io/EpiNow2/dev/reference/bound_dist.md)

- shape, scale:

  shape and scale parameters. Must be positive, `scale` strictly.

- rate:

  an alternative way to specify the scale.

- pmf:

  Probability mass of the given distribution; this is passed as a
  zero-indexed numeric vector (i.e. the fist entry represents the
  probability mass of zero). If not summing to one it will be normalised
  to sum to one internally.

- value:

  Value of the fixed (delta) distribution

## Value

A `dist_spec` representing a distribution of the given specification.

## Details

Probability distributions are ubiquitous in EpiNow2, usually
representing epidemiological delays (e.g., the generation time for
delays between becoming infecting and infecting others; or reporting
delays)

They are generated using functions that have a name corresponding to the
probability distribution that is being used. They generated `dist_spec`
objects that are then passed to the models underlying EpiNow2. All
parameters can be given either as fixed values (a numeric value) or as
uncertain values (a `dist_sepc`). If given as uncertain values,
currently only normally distributed parameters (generated using
`Normal()`) are supported.

Each distribution has a representation in terms of "natural" parameters
(the ones used in stan) but can sometimes also be specified using other
parameters such as the mean or standard deviation of the distribution.
If not given as natural parameters then these will be calculated from
the given parameters. If they have uncertainty, this will be done by
random sampling from the given uncertainty and converting resulting
parameters to their natural representation.

Currently available distributions are lognormal, gamma, normal, fixed
(delta) and nonparametric. The nonparametric is a special case where the
probability mass function is given directly as a numeric vector.

## See also

[`vignette("delays")`](https://epiforecasts.io/EpiNow2/dev/articles/delays.md)
for an overview of how delay distributions are used in EpiNow2

## Examples

``` r
LogNormal(mean = 4, sd = 1)
#> - lognormal distribution:
#>   meanlog:
#>     1.4
#>   sdlog:
#>     0.25
LogNormal(mean = 4, sd = 1, max = 10)
#> - lognormal distribution (max: 10):
#>   meanlog:
#>     1.4
#>   sdlog:
#>     0.25
# If specifying uncertain parameters, use the natural parameters
LogNormal(meanlog = Normal(1.5, 0.5), sdlog = 0.25, max = 10)
#> - lognormal distribution (max: 10):
#>   meanlog:
#>     - normal distribution:
#>       mean:
#>         1.5
#>       sd:
#>         0.5
#>   sdlog:
#>     0.25
Gamma(mean = 4, sd = 1)
#> - gamma distribution:
#>   shape:
#>     16
#>   rate:
#>     4
Gamma(shape = 16, rate = 4)
#> - gamma distribution:
#>   shape:
#>     16
#>   rate:
#>     4
Gamma(shape = Normal(16, 2), rate = Normal(4, 1))
#> - gamma distribution:
#>   shape:
#>     - normal distribution:
#>       mean:
#>         16
#>       sd:
#>         2
#>   rate:
#>     - normal distribution:
#>       mean:
#>         4
#>       sd:
#>         1
Normal(mean = 4, sd = 1)
#> - normal distribution:
#>   mean:
#>     4
#>   sd:
#>     1
Normal(mean = 4, sd = 1, max = 10)
#> - normal distribution (max: 10):
#>   mean:
#>     4
#>   sd:
#>     1
Exp(rate = 1)
#> - exp distribution:
#>   rate:
#>     1
Exp(mean = 4)
#> - exp distribution:
#>   rate:
#>     0.25
Weibull(shape = 1, scale = 1)
#> - weibull distribution:
#>   shape:
#>     1
#>   scale:
#>     1
Weibull(shape = 1, scale = 1, max = 10)
#> - weibull distribution (max: 10):
#>   shape:
#>     1
#>   scale:
#>     1
Weibull(mean = 4, sd = 1)
#> - weibull distribution:
#>   shape:
#>     4.5
#>   scale:
#>     4.4
NonParametric(c(0.1, 0.3, 0.2, 0.4))
#> - nonparametric distribution
#>   PMF: [0.1 0.3 0.2 0.4]
NonParametric(c(0.1, 0.3, 0.2, 0.1, 0.1))
#> - nonparametric distribution
#>   PMF: [0.12 0.37 0.25 0.12 0.12]
Fixed(value = 3)
#> - fixed value:
#>   3
Fixed(value = 3.5)
#> - fixed value:
#>   3.5
```
