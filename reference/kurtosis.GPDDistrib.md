# Excess Kurtosis of the Generalized Pareto Distribution

Closed form: \\\gamma_2 = 3(1-2\xi)(2\xi^2+\xi+3)/\\(1-3\xi)(1-4\xi)\\ -
3\\ for \\\xi \< 1/4\\, and `Inf` at or above one quarter. The threshold
is the tightest of the four moments, so this is the first quantity to
become unreportable as a fitted shape rises.

## Arguments

- x:

  A `GPDDistrib`, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value), each a numeric vector of length 1 or `n`. Only `xi` enters the
  value; settings with \\\xi \ge 1/4\\ give `Inf`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$sigma), length(theta$xi))`, `Inf` wherever \\\xi \ge
1/4\\.

## Notation

\\\xi\\ is the shape, of either sign. The scale does not enter a
standardized moment.

## See also

[`skewness.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GPDDistrib.md),
whose threshold is \\1/3\\;
[`kurtosis.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.ExponentialDistrib.md),
the \\\xi = 0\\ case;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## Examples

``` r
d <- gpd_distrib()

# Six at shape zero, and infinite at or above one quarter.
kurtosis(d, list(sigma = 1, xi = c(0, 0.2, 0.25, 0.3)))
#> [1]  6.0 70.8  Inf  Inf

# The four thresholds, from the mean's 1 down to this one's 1/4.
th <- list(sigma = 1, xi = 0.3)
c(mean = mean(d, th), variance = variance(d, th),
  skewness = skewness(d, th), kurtosis = kurtosis(d, th))
#>      mean  variance  skewness  kurtosis 
#>  1.428571  5.102041 16.443844       Inf 
```
