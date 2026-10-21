# Skewness of the Generalized Pareto Distribution

The skewness has the closed form \\\gamma_1 =
2(1+\xi)\sqrt{1-2\xi}/(1-3\xi)\\ for \\\xi \< 1/3\\, and is `Inf` at or
above one third. The scale cancels, so the value depends only on the
shape; it is 2 at \\\xi = 0\\, where the family is exponential, and it
climbs without bound as the shape approaches its threshold.

## Arguments

- x:

  A `GPDDistrib`, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value); each is a numeric vector of length 1 or `n`. Only `xi` enters
  the value; settings with \\\xi \ge 1/3\\ give `Inf`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$sigma), length(theta$xi))`, `Inf` wherever \\\xi \ge
1/3\\.

## Notation

\\\xi\\ is the shape, of either sign. The scale does not enter a
standardized moment.

## See also

[`kurtosis.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.GPDDistrib.md),
whose threshold is \\1/4\\;
[`skewness.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.ExponentialDistrib.md),
the \\\xi = 0\\ case;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

## Examples

``` r
d <- gpd_distrib()

# The skewness is two at shape zero, climbs steeply, and is infinite at or
# above one third.
skewness(d, list(sigma = 1, xi = c(0, 0.3, 1 / 3, 0.4)))
#> [1]  2.00000 16.44384      Inf      Inf
```
