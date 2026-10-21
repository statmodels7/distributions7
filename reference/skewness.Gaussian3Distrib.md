# Skewness of the Gaussian Distribution in Location and Precision

The skewness is exactly zero at every parameter value, because the
density is symmetric about \\\mu\\. A standardized moment is invariant
under a change of parametrization, so this is the same zero that the
other two Gaussian classes return.

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

[`kurtosis.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian3Distrib.md),
also zero;
[`skewness.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gaussian1Distrib.md);
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

## Examples

``` r
d <- gaussian3_distrib()

# The skewness is zero at every location and precision, with one value per
# setting.
skewness(d, list(mu = c(-2, 0, 5), tau = c(1, 4, 9)))
#> [1] 0 0 0
```
