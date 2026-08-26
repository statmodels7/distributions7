# Variance of the Gaussian Distribution

Closed form: \\\operatorname{Var}(Y) = \sigma^2\\. This parametrization
carries the standard deviation, so the variance is its square and
[`std_dev.distrib()`](https://statmodels7.github.io/distributions7/reference/std_dev.distrib.md)
returns the parameter back unchanged.

## Arguments

- x:

  A `Gaussian1Distrib`, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- theta:

  A named list with components `mu` (the mean) and `sigma` (the standard
  deviation, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`mean.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian1Distrib.md);
[`variance.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian2Distrib.md)
and
[`variance.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian3Distrib.md),
the same law written by its variance and by its precision;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

## Examples

``` r
d <- gaussian1_distrib()

# The square of the second parameter.
variance(d, list(mu = 0, sigma = c(1, 2, 4)))
#> [1]  1  4 16

# The standard deviation is the parameter itself.
all.equal(std_dev(d, list(mu = 0, sigma = 3)), 3)
#> [1] TRUE
```
