# Excess Kurtosis of the Negative Binomial Distribution in the NB1 Parametrization

Closed form: \\\gamma_2 = 6\theta/\mu + 1/\\\mu(1+\theta)\\\\, the
excess over the Gaussian. Both terms are positive, so the family is
always leptokurtic, and the expression tends to the Poisson's \\1/\mu\\
as the dispersion goes to zero.

## Arguments

- x:

  A `NegBin1Distrib`, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `theta` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$theta))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean and \\\theta \> 0\\ the dispersion.

## See also

[`skewness.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin1Distrib.md);
[`kurtosis.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin2Distrib.md),
the other parametrization;
[`kurtosis.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PoissonDistrib.md),
the limit;
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

## Examples

``` r
d <- negbin1_distrib()

# The published form, written out.
all.equal(kurtosis(d, list(mu = 4, theta = 2)), 12 / 4 + 1 / 12)
#> [1] TRUE

# It falls onto the Poisson's 1 / mu = 0.25 as the dispersion vanishes.
kurtosis(d, list(mu = 4, theta = c(1, 0.01, 1e-8)))
#> [1] 1.6250000 0.2625248 0.2500000
```
