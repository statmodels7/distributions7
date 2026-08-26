# Excess Kurtosis of the Inverse Gaussian Distribution

Closed form: \\\gamma_2 = 15\phi\mu\\, the excess over the Gaussian.
Like the skewness it depends on the two parameters through their product
alone, and it is exactly \\5\gamma_1^2/3\\, so the family occupies one
curve of the skewness-kurtosis plane, as the gamma does with its own
constant.

## Arguments

- x:

  An `InvGauss1Distrib`, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `phi` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$phi))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean and \\\phi \> 0\\ the dispersion.

## See also

[`skewness.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.InvGauss1Distrib.md),
to which this is tied;
[`kurtosis.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma2Distrib.md),
whose constant is \\3/2\\;
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

## Examples

``` r
d <- invgauss1_distrib()

# Fifteen times the product.
all.equal(kurtosis(d, list(mu = 2, phi = 0.5)), 15)
#> [1] TRUE

# The family lies on the curve kurtosis = (5/3) skewness^2.
th <- list(mu = 2, phi = 0.5)
all.equal(kurtosis(d, th), (5 / 3) * skewness(d, th)^2)
#> [1] TRUE
```
