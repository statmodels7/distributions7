# Skewness of the Laplace Distribution in Location and Rate

Exactly zero at every parameter value, the density being symmetric about
\\\mu\\. A standardized moment is free of the parametrization as well as
of the location and the scale, so this agrees with
[`skewness.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.LaplaceDistrib.md)
identically.

## Arguments

- x:

  A `Laplace2Distrib`, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or `n`. The values are not read, only their lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$lambda))`.

## See also

[`kurtosis.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Laplace2Distrib.md),
[`skewness.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.LaplaceDistrib.md),
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

## Examples

``` r
d <- laplace2_distrib()

# Zero at every location and rate, with one value per setting.
skewness(d, list(mu = c(0, 5), lambda = c(1, 2)))
#> [1] 0 0
```
