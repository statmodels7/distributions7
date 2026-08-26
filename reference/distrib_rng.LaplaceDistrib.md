# Laplace Random Number Generator

Draws `n` independent Laplace variates by **inverse transform**: `n`
uniform variates from
[`stats::runif()`](https://rdrr.io/r/stats/Uniform.html) are passed
through
[`distrib_quantile.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LaplaceDistrib.md).
The quantile function is elementary and exact, so this is both cheaper
and more accurate than the generalized ratio-of-uniforms fallback the
base class supplies, which would have to reject draws around the kink.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- n:

  A single positive integer, the number of draws. One uniform variate is
  consumed per draw.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled.
  `sigma` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_quantile.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LaplaceDistrib.md),
through which the draws are made;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back; and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()

# Inverse transform, so the draws are the quantiles of the uniforms drawn.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.5))
set.seed(2)
identical(a, distrib_quantile(d, runif(3), list(mu = 0.4, sigma = 1.5)))
#> [1] TRUE

# The sample variance recovers 2 sigma^2, not sigma^2.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
c(mean = mean(z), var = var(z), two_sigma_sq = 2 * 4)
#>         mean          var two_sigma_sq 
#>     3.007114     7.982085     8.000000 
```
