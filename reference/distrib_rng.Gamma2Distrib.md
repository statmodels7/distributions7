# Gamma Random Number Generator in Mean and Variance

Draws `n` independent gamma variates by calling
[`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html) at shape
\\\alpha = \mu^2/\sigma^2\\ and rate \\\lambda = \mu/\sigma^2\\, so the
draws come from R's own gamma generator and depend on `.Random.seed` in
the usual way. The ratio-of-uniforms fallback the base class supplies is
bypassed.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one variate per parameter setting. Both
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` strictly positive draws.

## See also

[`distrib_quantile.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gamma2Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()

# Same generator as stats::rgamma at the implied shape and rate.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 3, sigma2 = 2))
set.seed(2)
identical(a, rgamma(3, shape = 9 / 2, rate = 3 / 2))
#> [1] TRUE

# The sample moments recover both parameters directly.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 3, sigma2 = 2))
c(mu = mean(z), sigma2 = var(z))
#>       mu   sigma2 
#> 2.995309 1.985772 
```
