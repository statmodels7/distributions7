# Exponential Random Number Generator

Draws `n` independent exponential variates by calling
[`stats::rexp()`](https://rdrr.io/r/stats/Exponential.html) at
`rate = 1/mu`, so the draws come from R's own generator and depend on
`.Random.seed` in the usual way. The generalized ratio-of-uniforms
fallback the base class supplies is bypassed.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of length `n`. A value of length 1 is recycled, so a
  vector of length `n` draws one variate per mean. `mu` must be strictly
  positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` non-negative draws.

## See also

[`distrib_quantile.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ExponentialDistrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the mean back, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()

# Same generator as stats::rexp, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 2))
set.seed(2)
identical(a, rexp(3, rate = 1 / 2))
#> [1] TRUE

# The sample mean and standard deviation both estimate mu.
set.seed(21)
z <- distrib_rng(d, 2e4, list(mu = 3))
c(mean = mean(z), sd = sd(z))
#>     mean       sd 
#> 2.946014 2.970838 
```
