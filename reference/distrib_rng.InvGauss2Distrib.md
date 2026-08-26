# Inverse Gaussian Random Number Generator in Mean and Shape

Draws `n` independent inverse Gaussian variates by calling
[`statmod::rinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = 1/lambda`, so the draws come from that
package's generator and depend on `.Random.seed` in the usual way. The
ratio-of-uniforms fallback the base class supplies is bypassed.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one variate per parameter setting. Both
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` strictly positive draws.

## See also

[`distrib_quantile.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGauss2Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()

# Same generator as statmod::rinvgauss at dispersion 1/lambda.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 2, lambda = 3))
set.seed(2)
identical(a, statmod::rinvgauss(3, mean = 2, dispersion = 1 / 3))
#> [1] TRUE

# Both maximum likelihood estimates are closed form: the mean is the sample
# mean and the shape the reciprocal of mean(1/y) - 1/ybar.
set.seed(3)
z <- distrib_rng(d, 2e4, list(mu = 2, lambda = 3))
c(mu = mean(z), lambda = 1 / (mean(1 / z) - 1 / mean(z)))
#>       mu   lambda 
#> 1.994885 2.953600 
```
