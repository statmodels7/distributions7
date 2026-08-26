# Excess Kurtosis of the Negative Binomial Distribution

Closed form, replacing the numerical default: \\6/\theta +
\theta/\\\mu(\theta+\mu)\\\\, the excess over the Gaussian's three. Both
terms are positive, so the family is always leptokurtic, and the second
is what survives in the Poisson limit, where the expression tends to
\\1/\mu\\.

## Arguments

- x:

  A `NegBin2Distrib`, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `theta`
  (the dispersion, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$theta))`.

## See also

[`skewness.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin2Distrib.md),
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
[`kurtosis.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PoissonDistrib.md)
for the limit.

## Examples

``` r
d <- negbin2_distrib()

# The published form, written out.
mu <- 4; th <- 2
all.equal(kurtosis(d, list(mu = mu, theta = th)),
          6 / th + th / (mu * (th + mu)))
#> [1] TRUE

# The Poisson limit is 1 / mu = 0.25.
kurtosis(d, list(mu = 4, theta = c(1, 100, 1e8)))
#> [1] 6.0500000 0.3003846 0.2500001
```
