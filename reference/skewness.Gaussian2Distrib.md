# Skewness of the Gaussian Distribution in Location and Variance

Exactly zero at every parameter value, the density being symmetric about
\\\mu\\. A standardized moment does not see a change of parametrization,
so this is the same zero
[`skewness.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gaussian1Distrib.md)
returns.

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

[`kurtosis.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian2Distrib.md),
also zero;
[`skewness.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gaussian1Distrib.md);
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Examples

``` r
d <- gaussian2_distrib()

# Zero at every location and variance, with one value per setting.
skewness(d, list(mu = c(-2, 0, 5), sigma2 = c(1, 4, 9)))
#> [1] 0 0 0
```
