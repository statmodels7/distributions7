# Mean of the Negative Binomial Distribution in the NB1 Parametrization

Closed form: \\E\[Y\] = \mu\\. This parametrization carries the mean as
its first parameter, so the method reads it off. What distinguishes NB1
from NB2 is the variance function and not the mean.

## Arguments

- x:

  A `NegBin1Distrib`, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `theta`
  (the dispersion, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$theta))`.

## See also

[`variance.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin1Distrib.md),
which is linear in the mean;
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
which is quadratic;
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

## Examples

``` r
d <- negbin1_distrib()

# The first parameter is the mean, and the dispersion does not move it.
mean(d, list(mu = c(1, 2, 3), theta = 2))
#> [1] 1 2 3
```
