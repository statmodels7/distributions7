# Variance of the Laplace Distribution

The closed form that replaces the numerical default is
\\\operatorname{Var}(Y) = 2\sigma^2\\, so the standard deviation is
\\\sqrt2\\\sigma\\ and the scale is smaller than it by the factor
\\\sqrt2\\.

## Arguments

- x:

  A `LaplaceDistrib`, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- theta:

  A named list with components `mu` (the location) and `sigma` (the
  scale, positive); each is a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`. The location does not
enter the value, so a setting that varies only `mu` repeats one number.

## See also

[`mean.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.LaplaceDistrib.md),
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
[`variance.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Laplace2Distrib.md)
for the rate parametrization,
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()

# The variance is twice the square of the scale.
all.equal(variance(d, list(mu = 0, sigma = 2)), 8)
#> [1] TRUE

# The standard deviation is sqrt(2) times the scale.
std_dev(d, list(mu = 0, sigma = 1)) / 1
#> [1] 1.414214
```
