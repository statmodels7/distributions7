# Excess Kurtosis of the Laplace Distribution

Exactly 3 at every parameter value, the excess over the Gaussian's own
fourth standardized moment. A standardized moment is free of location
and scale, and the Laplace has no shape parameter, so the whole family
sits at one point of the kurtosis axis. That point is the reference the
toolkit's heavier-tailed families are read against.

## Arguments

- x:

  A `LaplaceDistrib`, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or `n`. The values are not read, only their lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of 3s, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`skewness.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.LaplaceDistrib.md),
[`variance.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.LaplaceDistrib.md),
[`kurtosis.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PseudoHuberDistrib.md),
which approaches this value as its shape goes to zero;
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()

# Three, whatever the location and the scale.
kurtosis(d, list(mu = c(0, 5), sigma = c(1, 5)))
#> [1] 3 3
```
