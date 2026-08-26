# Variance of the Inverse Gaussian Distribution

Closed form: \\\operatorname{Var}(Y) = \phi\mu^3\\. The cube is what
distinguishes this family from the gamma, whose variance function is
\\\mu^2\\: the inverse Gaussian's spread grows faster with the mean, so
it is the heavier-tailed of the two standard positive-response families.

## Arguments

- x:

  An `InvGauss1Distrib`, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `phi` (the
  dispersion, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$phi))`.

## Notation

\\\mu \> 0\\ is the mean and \\\phi \> 0\\ the dispersion.

## See also

[`mean.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.InvGauss1Distrib.md),
[`skewness.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.InvGauss1Distrib.md);
[`variance.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma2Distrib.md),
whose variance function is the square;
[`variance.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss2Distrib.md)
for the other parametrization;
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

## Examples

``` r
d <- invgauss1_distrib()

# The dispersion times the cube of the mean.
all.equal(variance(d, list(mu = 2, phi = 0.5)), 0.5 * 2^3)
#> [1] TRUE

# The variance function is cubic, so it climbs faster than a gamma's.
variance(d, list(mu = c(1, 2, 4), phi = 0.5))
#> [1]  0.5  4.0 32.0
```
