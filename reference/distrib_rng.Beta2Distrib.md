# Beta Random Number Generator in the Shapes

Draws `n` independent beta variates by calling
[`stats::rbeta()`](https://rdrr.io/r/stats/Beta.html) at
`shape1 = alpha` and `shape2 = beta`, so the draws come from R's own
beta generator and depend on `.Random.seed` in the usual way. The
ratio-of-uniforms fallback the base class supplies is bypassed.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one variate per parameter setting. Both
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws in \\(0, 1)\\.

## See also

[`distrib_quantile.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Beta2Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()

# Same generator as stats::rbeta at these two shapes.
set.seed(2)
a <- distrib_rng(d, 3, list(alpha = 2, beta = 5))
set.seed(2)
identical(a, rbeta(3, 2, 5))
#> [1] TRUE

# The sample moments recover the shapes through the concentration
# k = mean(1 - mean)/var - 1, with alpha = k mean and beta = k(1 - mean).
set.seed(8)
z <- distrib_rng(d, 2e4, list(alpha = 2, beta = 5))
m <- mean(z)
k <- m * (1 - m) / var(z) - 1
c(alpha = k * m, beta = k * (1 - m))
#>    alpha     beta 
#> 1.983183 4.951688 
```
