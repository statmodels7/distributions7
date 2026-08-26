# Mean of the Gamma Distribution

Closed form: \\E\[Y\] = \mu\\. This parametrization carries the mean and
the variance as its two parameters, so both of the first two moments are
reads, and the shape and rate are recovered from them where they are
needed.

## Arguments

- x:

  A `Gamma2Distrib`, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `sigma2`
  (the variance, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma2))`.

## See also

[`variance.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma2Distrib.md),
the other parameter;
[`skewness.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma2Distrib.md);
[`mean.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gamma1Distrib.md)
for the mean-dispersion parametrization;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

## Examples

``` r
d <- gamma2_distrib()

# The first parameter is the mean.
mean(d, list(mu = c(1, 2, 3), sigma2 = 1))
#> [1] 1 2 3
```
