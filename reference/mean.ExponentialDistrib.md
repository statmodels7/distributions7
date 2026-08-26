# Mean of the Exponential Distribution

Closed form: \\E\[Y\] = \mu\\. The family carries one parameter and it
is the mean, so this method reads it off. The rate is its reciprocal.

## Arguments

- x:

  An `ExponentialDistrib`, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- theta:

  A named list with one component, `mu` (the mean, positive), a numeric
  vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, the length of `theta$mu`.

## See also

[`variance.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.ExponentialDistrib.md),
which is the square of this;
[`mean.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Weibull1Distrib.md)
and
[`mean.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gamma2Distrib.md),
the two families that contain this one;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

## Examples

``` r
d <- exponential_distrib()

# The one parameter is the mean.
mean(d, list(mu = c(1, 2, 3)))
#> [1] 1 2 3
```
