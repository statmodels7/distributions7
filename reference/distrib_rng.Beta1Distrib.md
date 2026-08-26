# Beta Random Number Generator in Mean and Precision

Draws `n` independent beta variates by calling
[`stats::rbeta()`](https://rdrr.io/r/stats/Beta.html) at shapes \\\alpha
= \mu\phi\\ and \\\beta = (1-\mu)\phi\\, so the draws come from R's own
beta generator and depend on `.Random.seed` in the usual way. The
ratio-of-uniforms fallback the base class supplies is bypassed.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of length `n`. A component of length 1 is recycled, so a
  vector of length `n` draws one variate per parameter setting. `mu`
  must lie strictly in \\(0, 1)\\ and `phi` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws in \\(0, 1)\\.

## See also

[`distrib_quantile.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Beta1Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()

# Same generator as stats::rbeta at the implied shapes.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0.4, phi = 5))
set.seed(2)
identical(a, rbeta(3, 2, 3))
#> [1] TRUE

# The sample moments recover the parameters: the mean directly, and the
# precision as mu(1 - mu)/var - 1.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 0.4, phi = 5))
c(mu = mean(z), phi = mean(z) * (1 - mean(z)) / var(z) - 1)
#>        mu       phi 
#> 0.4000021 5.0221098 
```
