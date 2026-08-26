# Mean of the Gaussian Distribution in Location and Precision

Closed form: \\E\[Y\] = \mu\\. This is the Gaussian written by its
precision \\\tau = 1/\sigma^2\\, the parametrization a Bayesian prior or
a penalized likelihood usually reaches for, and the mean is the location
as it is in the other two.

## Arguments

- x:

  A `Gaussian3Distrib`, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- theta:

  A named list with components `mu` (the mean, any real value) and `tau`
  (the precision, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$tau))`.

## See also

[`variance.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian3Distrib.md),
the reciprocal of the precision;
[`mean.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian1Distrib.md)
and
[`mean.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian2Distrib.md);
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

## Examples

``` r
d <- gaussian3_distrib()

# The first parameter is the mean, and the precision does not move it.
mean(d, list(mu = c(-1, 0, 4), tau = 2))
#> [1] -1  0  4
```
