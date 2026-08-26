# Excess Kurtosis of the Gumbel Distribution

Constant: \\\gamma_2 = 12/5 = 2.4\\, the excess over the Gaussian. Like
the skewness it is fixed for the whole family, there being no shape
parameter for a standardized moment to depend on.

## Arguments

- x:

  A `GumbelDistrib`, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or `n`. The values are not read, only their lengths.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of 2.4s, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`skewness.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GumbelDistrib.md),
the other constant;
[`variance.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.GumbelDistrib.md);
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

## Examples

``` r
d <- gumbel_distrib()

# 12 / 5 for every Gumbel.
kurtosis(d, list(mu = c(0, 5), sigma = c(1, 7)))
#> [1] 2.4 2.4
```
