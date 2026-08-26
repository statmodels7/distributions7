# Variance of the Gamma Distribution in Mean and Dispersion

Closed form: \\\operatorname{Var}(Y) = \phi\mu^2\\. The quadratic
variance function is why the gamma is the standard family for a positive
response whose spread grows in proportion to its level: the coefficient
of variation is \\\sqrt\phi\\, the same at every mean.

## Arguments

- x:

  A `Gamma1Distrib`, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `phi` (positive),
  each a numeric vector of length 1 or `n`. Small `phi` is a large
  shape, where the family approaches a Gaussian.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$phi))`.

## Notation

\\\mu \> 0\\ is the mean and \\\phi \> 0\\ the dispersion. The implied
shape is \\1/\phi\\.

## See also

[`mean.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gamma1Distrib.md),
[`skewness.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma1Distrib.md);
[`variance.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma2Distrib.md),
the same quantity as a parameter;
[`variance.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss1Distrib.md),
whose variance function is cubic;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

## Examples

``` r
d <- gamma1_distrib()

# The dispersion times the square of the mean.
all.equal(variance(d, list(mu = 2, phi = 0.25)), 1)
#> [1] TRUE

# The coefficient of variation is sqrt(phi) at every mean.
th <- list(mu = c(1, 10, 100), phi = 0.25)
std_dev(d, th) / mean(d, th)
#> [1] 0.5 0.5 0.5
```
