# Excess Kurtosis of the Gaussian Distribution in Location and Variance

Exactly zero at every parameter value.
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
reports the excess over the Gaussian's own fourth standardized moment of
3, so every Gaussian sits at the origin of the scale whichever of the
three parametrizations it is written in.

## Arguments

- x:

  A `Gaussian2Distrib`, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- theta:

  A named list with components `mu` and `sigma2` (positive), each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$sigma2))`.

## See also

[`skewness.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gaussian2Distrib.md),
also zero;
[`kurtosis.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian1Distrib.md);
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Examples

``` r
d <- gaussian2_distrib()

# Zero for every Gaussian, in every parametrization.
kurtosis(d, list(mu = 0, sigma2 = c(1, 100)))
#> [1] 0 0
```
