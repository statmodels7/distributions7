# Logistic Random Number Generator

Draws `n` independent logistic variates by calling
[`stats::rlogis()`](https://rdrr.io/r/stats/Logistic.html), so the draws
come from R's own generator and depend on `.Random.seed` in the usual
way. The generalized ratio-of-uniforms fallback the base class supplies
is bypassed.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled.
  `sigma` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_quantile.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LogisticDistrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()

# Same generator as stats::rlogis, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.5))
set.seed(2)
identical(a, rlogis(3, location = 0.4, scale = 1.5))
#> [1] TRUE

# The sample variance recovers pi^2 sigma^2 / 3, not sigma^2.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
c(mean = mean(z), var = var(z), pi_sq_over_3 = pi^2 * 4 / 3)
#>         mean          var pi_sq_over_3 
#>     3.009384    13.135972    13.159473 
```
