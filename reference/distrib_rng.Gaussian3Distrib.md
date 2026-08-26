# Gaussian Random Number Generator in Mean and Precision

Draws `n` independent Gaussian variates by calling
[`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html) at
`sd = 1/sqrt(tau)`, so the draws come from R's own normal generator and
depend on `.Random.seed` in the usual way. The inverse-transform
fallback the base class supplies is bypassed.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of length `n`. A component of length 1 is recycled, so a
  vector of length `n` draws one variate per parameter setting. `tau`
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_quantile.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gaussian3Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()

# Same generator as stats::rnorm at sd = 1/sqrt(tau).
set.seed(2)
a <- distrib_rng(d, 3, list(mu = c(0, 1, 2), tau = c(1, 0.25, 0.0625)))
set.seed(2)
identical(a, rnorm(3, mean = c(0, 1, 2), sd = c(1, 2, 4)))
#> [1] TRUE

# The sample moments recover the parameters, the precision as 1/var.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 3, tau = 0.25))
c(mean = mean(z), tau = 1 / var(z))
#>      mean       tau 
#> 3.0104856 0.2471227 
```
