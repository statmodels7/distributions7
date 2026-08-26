# Mean of the Laplace Distribution

Closed form, replacing the numerical default: \\E\[Y\] = \mu\\. The
density is symmetric about \\\mu\\, so the location is the mean and the
median at once. The value is recycled to the length the two parameters
imply, and no quadrature is run: the kink at \\y = \mu\\ makes the
numerical route needlessly awkward for a quantity available in one read.

## Arguments

- x:

  A `LaplaceDistrib`, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- theta:

  A named list with components `mu` (the location) and `sigma` (the
  scale, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`variance.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.LaplaceDistrib.md),
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
[`mean.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Laplace2Distrib.md)
for the same family written by its rate,
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()

# The location is the mean, and the scale does not move it.
mean(d, list(mu = c(-1, 0, 4), sigma = 3))
#> [1] -1  0  4
```
