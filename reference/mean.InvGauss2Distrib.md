# Mean of the Inverse Gaussian Distribution in Mean and Rate

Closed form: \\E\[Y\] = \mu\\. This is the inverse Gaussian written by
its mean and the rate \\\lambda = 1/\phi\\, the parametrization the
classical literature uses, and the mean is the first parameter in
either.

## Arguments

- x:

  An `InvGauss2Distrib`, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `lambda`
  (the rate, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$lambda))`.

## See also

[`variance.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss2Distrib.md),
where the rate does enter;
[`mean.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.InvGauss1Distrib.md)
for the dispersion parametrization;
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

## Examples

``` r
d <- invgauss2_distrib()

# The first parameter is the mean, and the rate does not move it.
mean(d, list(mu = c(1, 2, 3), lambda = 8))
#> [1] 1 2 3
```
