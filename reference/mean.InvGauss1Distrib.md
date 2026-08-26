# Mean of the Inverse Gaussian Distribution

Closed form: \\E\[Y\] = \mu\\. This parametrization carries the mean and
a dispersion, so the mean is a read and the higher moments are functions
of the product \\\phi\mu\\.

## Arguments

- x:

  An `InvGauss1Distrib`, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `phi` (the
  dispersion, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$phi))`.

## See also

[`variance.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss1Distrib.md),
which is \\\phi\mu^3\\;
[`skewness.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.InvGauss1Distrib.md);
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

## Examples

``` r
d <- invgauss1_distrib()

# The first parameter is the mean.
mean(d, list(mu = c(1, 2, 3), phi = 0.5))
#> [1] 1 2 3
```
