# Skewness of the Negative Binomial Distribution

Closed form, replacing the numerical default: \\(\theta +
2\mu)/\sqrt{\mu\theta(\theta+\mu)}\\. It is positive at every parameter
value, a count distribution having a right tail and a floor at zero, and
it decreases towards the Poisson value \\\mu^{-1/2}\\ as \\\theta\\
grows.

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

A numeric vector, of length
`max(length(theta$mu), length(theta$theta))`.

## See also

[`kurtosis.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin2Distrib.md),
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
[`skewness.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.PoissonDistrib.md)
for the limit.

## Examples

``` r
d <- negbin2_distrib()

# The published form, written out.
mu <- 4; th <- 2
all.equal(skewness(d, list(mu = mu, theta = th)),
          (th + 2 * mu) / sqrt(mu * th * (th + mu)))
#> [1] TRUE

# It falls onto the Poisson's mu^(-1/2) = 0.5 as the dispersion grows.
skewness(d, list(mu = 4, theta = c(1, 100, 1e8)))
#> [1] 2.0124612 0.5295136 0.5000000
```
