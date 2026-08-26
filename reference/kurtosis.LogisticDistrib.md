# Excess Kurtosis of the Logistic Distribution

Constant: \\\gamma_2 = 6/5 = 1.2\\, the excess over the Gaussian. The
family is location-scale with no shape parameter, so a standardized
moment is a number: every logistic has the same tail weight, a little
heavier than a Gaussian's and well short of a Laplace's 3.

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

A numeric vector of 1.2s, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`skewness.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.LogisticDistrib.md),
the other constant;
[`kurtosis.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian1Distrib.md)
and
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
the two families it sits between;
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

## Examples

``` r
d <- logistic_distrib()

# 6 / 5 for every logistic.
kurtosis(d, list(mu = c(0, 5), sigma = c(1, 3)))
#> [1] 1.2 1.2

# Between the Gaussian's 0 and the Laplace's 3.
c(gaussian = kurtosis(gaussian1_distrib(), list(mu = 0, sigma = 1)),
  logistic = kurtosis(d, list(mu = 0, sigma = 1)),
  laplace  = kurtosis(laplace_distrib(), list(mu = 0, sigma = 1)))
#> gaussian logistic  laplace 
#>      0.0      1.2      3.0 
```
