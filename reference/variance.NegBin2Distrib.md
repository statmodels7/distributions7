# Variance of the Negative Binomial Distribution

Closed form, replacing the numerical default: \\\operatorname{Var}(Y) =
\mu + \mu^2/\theta\\. The quadratic term is the overdispersion, so the
variance exceeds the mean at every finite \\\theta\\ and falls back onto
it as \\\theta\\ grows, which is the Poisson limit of the family.

## Arguments

- x:

  A `NegBin2Distrib`, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `theta`
  (the dispersion, positive), each a numeric vector of length 1 or `n`.
  Small `theta` means heavy overdispersion; the variance diverges as it
  approaches zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$theta))`.

## See also

[`mean.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.NegBin2Distrib.md),
[`skewness.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin2Distrib.md),
[`kurtosis.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin2Distrib.md);
[`variance.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PoissonDistrib.md)
for the equidispersed limit.

## Examples

``` r
d <- negbin2_distrib()

# Mean 4, dispersion 2: the variance is 4 + 16/2.
all.equal(variance(d, list(mu = 4, theta = 2)), 12)
#> [1] TRUE

# The overdispersion vanishes as theta grows, leaving the Poisson variance.
variance(d, list(mu = 4, theta = c(1, 10, 1e4)))
#> [1] 20.0000  5.6000  4.0016
```
