# Mean of the Gamma Distribution in Mean and Dispersion

Closed form: \\E\[Y\] = \mu\\. This is the gamma written by its mean and
a dispersion, the parametrization a generalized linear model uses, so
the mean is the first parameter and the variance function is
\\\phi\mu^2\\.

## Arguments

- x:

  A `Gamma1Distrib`, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `phi` (the
  dispersion, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$phi))`.

## See also

[`variance.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma1Distrib.md),
which is \\\phi\mu^2\\;
[`mean.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gamma2Distrib.md)
for the mean-variance parametrization;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

## Examples

``` r
d <- gamma1_distrib()

# The first parameter is the mean, and the dispersion does not move it.
mean(d, list(mu = c(1, 2, 3), phi = 0.25))
#> [1] 1 2 3
```
