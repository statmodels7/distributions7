# Cauchy Random Number Generator

Draws `n` independent Cauchy variates by calling
[`stats::rcauchy()`](https://rdrr.io/r/stats/Cauchy.html), so the draws
come from R's own generator and depend on `.Random.seed` in the usual
way. The generalized ratio-of-uniforms fallback the base class supplies
is bypassed.

A sample from this family carries extreme values at a rate that does not
diminish with `n`, and its running mean does not settle: the sample mean
of `n` Cauchy draws is itself Cauchy with the same scale, whatever `n`
is.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

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

[`distrib_quantile.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.CauchyDistrib.md)
for the inverse-transform route,
[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md)
for the non-convergence of the sample mean, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()

# Same generator as stats::rcauchy, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.5))
set.seed(2)
identical(a, rcauchy(3, location = 0.4, scale = 1.5))
#> [1] TRUE

# The sample median converges on mu; the sample mean does not converge at
# all, being itself Cauchy at every sample size.
set.seed(4)
z <- distrib_rng(d, 1e5, list(mu = 3, sigma = 2))
n <- c(1e2, 1e3, 1e4, 1e5)
rbind(median = vapply(n, function(k) median(z[1:k]), numeric(1)),
      mean   = vapply(n, function(k) mean(z[1:k]), numeric(1)))
#>            [,1]     [,2]      [,3]      [,4]
#> median 2.693216 3.143942  3.000707  2.998429
#> mean   2.097306 4.772852 -1.155020 -3.675191
```
