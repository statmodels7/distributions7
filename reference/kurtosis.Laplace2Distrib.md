# Excess Kurtosis of the Laplace Distribution in Location and Rate

Exactly 3 at every parameter value, the excess over the Gaussian. A
standardized moment does not see a change of parametrization, so this is
the same number
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md)
returns and it moves with neither the location nor the rate.

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

A numeric vector of 3s, of length
`max(length(theta$mu), length(theta$lambda))`.

## See also

[`skewness.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Laplace2Distrib.md),
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

## Examples

``` r
d <- laplace2_distrib()

# Three, in either parametrization.
c(kurtosis(d, list(mu = 0, lambda = 2)),
  kurtosis(laplace_distrib(), list(mu = 0, sigma = 0.5)))
#> [1] 3 3
```
