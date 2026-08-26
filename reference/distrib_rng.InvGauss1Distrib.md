# Inverse Gaussian Random Number Generator in Mean and Dispersion

Draws `n` independent inverse Gaussian variates by calling
[`statmod::rinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = phi`, so the draws come from that
package's generator and depend on `.Random.seed` in the usual way. The
ratio-of-uniforms fallback the base class supplies is bypassed.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of length `n`. A component of length 1 is recycled, so a
  vector of length `n` draws one variate per parameter setting. Both
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` strictly positive draws.

## See also

[`distrib_quantile.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGauss1Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()

# Same generator as statmod::rinvgauss at this parametrization.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 1, phi = 2))
set.seed(2)
identical(a, statmod::rinvgauss(3, mean = 1, dispersion = 2))
#> [1] TRUE

# Both maximum likelihood estimates are closed form here: the mean is the
# sample mean and the dispersion the mean of 1/y minus 1/ybar.
set.seed(3)
z <- distrib_rng(d, 2e4, list(mu = 1, phi = 2))
c(mu = mean(z), phi = mean(1 / z) - 1 / mean(z))
#>       mu      phi 
#> 1.001860 2.031436 
```
