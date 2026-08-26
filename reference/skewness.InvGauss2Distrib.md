# Skewness of the Inverse Gaussian Distribution in Mean and Rate

Closed form: \\\gamma_1 = 3\sqrt{\mu/\lambda}\\. It is positive at every
parameter value and depends on the two only through their ratio, so the
whole family is one curve indexed by \\\mu/\lambda\\. Small mean or
large rate gives a nearly Gaussian shape.

## Arguments

- x:

  An `InvGauss2Distrib`, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- theta:

  A named list with components `mu` (positive) and `lambda` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$mu), length(theta$lambda))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean and \\\lambda \> 0\\ the rate.

## See also

[`kurtosis.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.InvGauss2Distrib.md),
which is \\5/3\\ times the square of this;
[`skewness.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.InvGauss1Distrib.md);
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

## Examples

``` r
d <- invgauss2_distrib()

# Three times the square root of the ratio.
all.equal(skewness(d, list(mu = 2, lambda = 8)), 1.5)
#> [1] TRUE

# It agrees with the dispersion parametrization on the same law.
all.equal(skewness(d, list(mu = 2, lambda = 8)),
          skewness(invgauss1_distrib(), list(mu = 2, phi = 1 / 8)))
#> [1] TRUE
```
