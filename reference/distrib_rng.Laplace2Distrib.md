# Laplace Random Number Generator, Rate Parametrization

Draws `n` independent Laplace variates by **inverse transform**: `n`
uniform variates from
[`stats::runif()`](https://rdrr.io/r/stats/Uniform.html) are passed
through
[`distrib_quantile.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Laplace2Distrib.md).
The quantile function is elementary and exact, so this is cheaper and
more accurate than the generalized ratio-of-uniforms fallback the base
class supplies.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- n:

  A single positive integer, the number of draws. One uniform variate is
  consumed per draw.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled.
  `lambda` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_quantile.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Laplace2Distrib.md),
through which the draws are made;
[`distrib_rng.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.LaplaceDistrib.md)
for the scale parametrization; and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()

# Inverse transform, so the draws are the quantiles of the uniforms drawn.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0.4, lambda = 2))
set.seed(2)
identical(a, distrib_quantile(d, runif(3), list(mu = 0.4, lambda = 2)))
#> [1] TRUE

# The sample variance recovers 2/lambda^2.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 3, lambda = 0.5))
c(mean = mean(z), var = var(z), two_over_lambda_sq = 2 / 0.5^2)
#>               mean                var two_over_lambda_sq 
#>           3.007114           7.982085           8.000000 
```
