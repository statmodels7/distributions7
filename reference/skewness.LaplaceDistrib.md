# Skewness of the Laplace Distribution

Exactly zero at every parameter value. The density is symmetric about
\\\mu\\, so every odd central moment vanishes. The constant is returned
directly, recycled to the length the parameters imply.

## Arguments

- x:

  A `LaplaceDistrib`, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or `n`. The values are not read, only their lengths, and
  the list is still aligned and validated.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
which is 3;
[`mean.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.LaplaceDistrib.md);
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()

# Zero at every location and scale, with one value per setting.
skewness(d, list(mu = c(0, 5), sigma = c(1, 2)))
#> [1] 0 0
```
