# Variance of the Poisson-Inverse Gaussian Distribution

Closed form, the second cumulant: \\\operatorname{Var}(Y) = \kappa_2 =
\mu + \sigma\mu^2\\. The quadratic term is the overdispersion the
inverse Gaussian mixing adds, and it has the same shape as a negative
binomial's; the two families separate at the third cumulant and above,
where this one is the heavier-tailed.

## Arguments

- x:

  A `Pig1Distrib`, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `sigma` (positive),
  each a numeric vector of length 1 or `n`. The Poisson limit is
  \\\sigma \to 0\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu \> 0\\ is the mean, \\\sigma \> 0\\ the dispersion and
\\\kappa_k\\ the \\k\\-th cumulant.

## See also

[`mean.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Pig1Distrib.md),
[`skewness.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig1Distrib.md);
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
which matches this shape;
[`variance.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PoissonDistrib.md),
the limit;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

## Examples

``` r
d <- pig1_distrib()

# Mean plus dispersion times the square of the mean.
all.equal(variance(d, list(mu = 3, sigma = 1)), 12)
#> [1] TRUE

# It falls onto the Poisson variance as the dispersion vanishes.
c(pig = variance(d, list(mu = 3, sigma = 1e-9)),
  poisson = variance(poisson_distrib(), list(mu = 3)))
#>     pig poisson 
#>       3       3 
```
