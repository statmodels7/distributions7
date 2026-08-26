# Skewness of the Negative Binomial Distribution in the NB1 Parametrization

Closed form: \\\gamma_1 = (1+2\theta)/\sqrt{\mu(1+\theta)}\\. It is
positive at every parameter value and tends to the Poisson's
\\\mu^{-1/2}\\ as the dispersion goes to zero, which is where NB1
approaches a Poisson.

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

A numeric vector, of length
`max(length(theta$mu), length(theta$theta))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean and \\\theta \> 0\\ the dispersion.

## See also

[`kurtosis.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin1Distrib.md);
[`skewness.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin2Distrib.md),
the other parametrization;
[`skewness.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.PoissonDistrib.md),
the limit;
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

## Examples

``` r
d <- negbin1_distrib()

# The published form, written out.
all.equal(skewness(d, list(mu = 4, theta = 2)), 5 / sqrt(4 * 3))
#> [1] TRUE

# It falls onto the Poisson's mu^(-1/2) = 0.5 as the dispersion vanishes.
skewness(d, list(mu = 4, theta = c(1, 0.01, 1e-8)))
#> [1] 1.060660 0.507469 0.500000
```
