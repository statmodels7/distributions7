# Student t Random Number Generator

Draws `n` independent location-scale Student t variates as \\\mu +
\sigma T\\ with \\T\\ from
[`stats::rt()`](https://rdrr.io/r/stats/TDist.html), so the draws come
from R's own generator and depend on `.Random.seed` in the usual way.
The generalized ratio-of-uniforms fallback the base class supplies is
bypassed.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of length `n`. A component of length 1 is
  recycled, so a vector of length `n` draws one variate per parameter
  setting. `sigma` and `nu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_quantile.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.StudentT1Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()

# The draws are stats::rt shifted and scaled, so the same seed reproduces
# them.
set.seed(12)
a <- distrib_rng(d, 3, list(mu = 0.4, sigma = 1.2, nu = 5))
set.seed(12)
identical(a, 0.4 + 1.2 * rt(3, df = 5))
#> [1] TRUE

# The sample variance recovers sigma^2 * nu / (nu - 2), not sigma^2.
set.seed(13)
z <- distrib_rng(d, 5e4, list(mu = 0.4, sigma = 1.2, nu = 5))
c(sample = var(z), theoretical = 1.2^2 * 5 / 3, sigma_sq = 1.2^2)
#>      sample theoretical    sigma_sq 
#>    2.416863    2.400000    1.440000 
```
