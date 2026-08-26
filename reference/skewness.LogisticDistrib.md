# Skewness of the Logistic Distribution

Exactly zero at every parameter value. The density is symmetric about
\\\mu\\, so every odd central moment vanishes. The constant is recycled
to the length the parameters imply.

## Arguments

- x:

  A `LogisticDistrib`, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of zeros, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`kurtosis.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LogisticDistrib.md),
which is 6/5;
[`mean.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.LogisticDistrib.md);
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

## Examples

``` r
d <- logistic_distrib()

# Zero at every location and scale, with one value per setting.
skewness(d, list(mu = c(0, 5), sigma = c(1, 2)))
#> [1] 0 0
```
