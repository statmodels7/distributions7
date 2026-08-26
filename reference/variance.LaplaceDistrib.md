# Variance of the Laplace Distribution

Closed form, replacing the numerical default: \\\operatorname{Var}(Y) =
2\sigma^2\\. The scale is not the standard deviation here: it is smaller
by \\\sqrt2\\, which is worth knowing before a Laplace scale is compared
with a Gaussian one.

## Arguments

- x:

  A `LaplaceDistrib`, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- theta:

  A named list with components `mu` (the location) and `sigma` (the
  scale, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`. The location does not
enter the value, so a setting that varies `mu` alone repeats one number.

## See also

[`mean.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.LaplaceDistrib.md),
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md),
[`variance.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Laplace2Distrib.md)
for the rate parametrization,
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()

# Twice the square of the scale.
all.equal(variance(d, list(mu = 0, sigma = 2)), 8)
#> [1] TRUE

# The standard deviation is sqrt(2) times the scale.
std_dev(d, list(mu = 0, sigma = 1)) / 1
#> [1] 1.414214
```
