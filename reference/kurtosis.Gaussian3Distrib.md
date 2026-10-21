# Excess Kurtosis of the Gaussian Distribution in Location and Precision

The excess kurtosis is exactly zero at every parameter value, because
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
reports the excess over the Gaussian's own fourth standardized moment.
The parametrization does not enter a standardized moment.

## Arguments

- x:

  A `Gaussian3Distrib`, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- theta:

  A named list with components `mu` and `tau` (positive); each is a
  numeric vector of length 1 or `n`. The values are checked against the
  parameter bounds, and only their lengths enter the result.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$tau))`.

## See also

[`skewness.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gaussian3Distrib.md),
also zero;
[`kurtosis.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian1Distrib.md);
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

## Examples

``` r
d <- gaussian3_distrib()

# The excess kurtosis is zero for every Gaussian, in every parametrization.
kurtosis(d, list(mu = 0, tau = c(0.01, 100)))
#> [1] 0 0
```
