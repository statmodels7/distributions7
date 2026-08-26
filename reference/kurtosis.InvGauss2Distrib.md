# Excess Kurtosis of the Inverse Gaussian Distribution in Mean and Rate

Closed form: \\\gamma_2 = 15\mu/\lambda\\, the excess over the Gaussian.
Like the skewness it depends on the two parameters through their ratio
alone, and it is exactly \\5\gamma_1^2/3\\, so the family occupies one
curve of the skewness-kurtosis plane, as it does in the other
parametrization.

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

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$lambda))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean and \\\lambda \> 0\\ the rate.

## See also

[`skewness.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.InvGauss2Distrib.md),
to which this is tied;
[`kurtosis.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.InvGauss1Distrib.md);
[`kurtosis.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma1Distrib.md),
whose constant is \\3/2\\;
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

## Examples

``` r
d <- invgauss2_distrib()

# Fifteen times the ratio.
all.equal(kurtosis(d, list(mu = 2, lambda = 8)), 3.75)
#> [1] TRUE

# The family lies on the curve kurtosis = (5/3) skewness^2.
th <- list(mu = 2, lambda = 8)
all.equal(kurtosis(d, th), (5 / 3) * skewness(d, th)^2)
#> [1] TRUE
```
