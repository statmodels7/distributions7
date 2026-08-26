# Skewness of the Gaussian Distribution

Exactly zero at every parameter value. The density is symmetric about
\\\mu\\, so every odd central moment vanishes. The constant is recycled
to the length the parameters imply, and no quadrature is run: the
numerical default returns a small non-zero number here, and an exact
zero is available.

## Arguments

- x:

  A `Gaussian1Distrib`, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths, and the list is still aligned and validated.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`kurtosis.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian1Distrib.md),
zero by the excess convention;
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the generic;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

## Examples

``` r
d <- gaussian1_distrib()

# Zero at every location and scale, with one value per setting.
skewness(d, list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#> [1] 0 0 0
```
