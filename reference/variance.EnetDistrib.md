# Variance of the Elastic-Net Distribution

Returns \\(1 + xG)/c\\ in closed form, with \\c = \lambda(1-\alpha)\\,
\\x = a/\sqrt c\\ and \\G = \mathrm{d}\log M/\mathrm{d}x\\. It is
\\-2\\\partial\log Z/\partial c\\, the first cumulant of the sufficient
statistic \\-z^2/2\\, so it needs nothing the score does not already
compute.

At \\\alpha \to 0\\ it tends to \\1/c\\, the Gaussian's; at \\\alpha \to
1\\ to \\2/a^2\\, the Laplace's.

## Arguments

- x:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- theta:

  A named list with components `mu`, `lambda` and `alpha`. Aligned and
  validated by name, so a missing or out-of-bounds component throws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, strictly positive, of the length the
recycled parameters imply. The location does not enter the value, so a
setting that varies it repeats one number.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, \\x = a/\sqrt c\\,
\\M\\ the Mills ratio and \\Z\\ the normalizing constant.

## See also

[`mean.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.EnetDistrib.md)
and
[`skewness.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.EnetDistrib.md),
and
[`distrib_expected_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.EnetDistrib.md),
which uses the same \\\log Z\\ derivatives.

## Examples

``` r
d <- enet_distrib()
th <- list(mu = 0, lambda = 2, alpha = 0.5)

# Against a quadrature of the second central moment.
c(closed = variance(d, th),
  quadrature = integrate(function(u) u^2 * distrib_pdf(d, u, th),
                         -Inf, Inf)$value)
#>     closed quadrature 
#>  0.4748647  0.4748647 

# The two ends: 1/c at alpha -> 0 and 2/a^2 at alpha -> 1.
rbind(alpha = c(1e-10, 1 - 1e-10),
      ours = c(variance(d, list(mu = 0, lambda = 2, alpha = 1e-10)),
               variance(d, list(mu = 0, lambda = 2, alpha = 1 - 1e-10))),
      limit = c(1 / 2, 2 / 2^2))
#>        [,1]      [,2]
#> alpha 1e-10 1.0000000
#> ours  5e-01 0.4999994
#> limit 5e-01 0.5000000
```
