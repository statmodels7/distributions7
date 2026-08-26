# Weibull Random Number Generator

Draws `n` independent Weibull variates by calling
[`stats::rweibull()`](https://rdrr.io/r/stats/Weibull.html), which
inverts the distribution function at a uniform draw. The inversion is
exact because the quantile function is elementary, so the generalized
ratio-of-uniforms fallback the base class supplies is bypassed. The
draws depend on `.Random.seed` in the usual way.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one variate per parameter setting. Both
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` positive draws.

## See also

[`distrib_quantile.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Weibull1Distrib.md)
for the function inverted,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()

# Same generator as stats::rweibull, so the same seed gives the same draws.
set.seed(3)
a <- distrib_rng(d, 3, list(mu = 2, sigma = 1.5))
set.seed(3)
identical(a, rweibull(3, shape = 1.5, scale = 2))
#> [1] TRUE

# The sample mean recovers mu * gamma(1 + 1/sigma), not mu itself.
set.seed(11)
z <- distrib_rng(d, 2e4, list(mu = 2, sigma = 1.5))
c(sample = mean(z), theoretical = 2 * gamma(1 + 1 / 1.5), scale = 2)
#>      sample theoretical       scale 
#>    1.792680    1.805491    2.000000 
```
