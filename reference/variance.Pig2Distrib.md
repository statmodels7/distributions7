# Variance of the Orthogonal Poisson-Inverse Gaussian Distribution

The second cumulant \\\kappa_2 = \mu + \sigma\mu^2\\ of
[`variance.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Pig1Distrib.md),
read at the dispersion \\\sigma = \\
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)`(mu, alpha)`
this parametrization implies. The law is the same; what changes is which
number is reported, and the orthogonal dispersion is the one whose
estimate is asymptotically uncorrelated with the mean's.

## Arguments

- x:

  A `Pig2Distrib`, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- theta:

  A named list with components `mu` (positive) and `alpha` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$alpha))`.

## Notation

\\\mu \> 0\\ is the mean, \\\alpha \> 0\\ the orthogonal dispersion,
\\\sigma\\ the dispersion of
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
and \\\kappa_k\\ the \\k\\-th cumulant.

## See also

[`mean.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Pig2Distrib.md),
[`skewness.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig2Distrib.md);
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the mapping;
[`variance.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Pig1Distrib.md);
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

## Examples

``` r
d <- pig2_distrib()

# The two parametrizations agree once the mapping is applied.
s <- distributions7:::pig2_sigma(3, 1)
all.equal(variance(d, list(mu = 3, alpha = 1)),
          variance(pig1_distrib(), list(mu = 3, sigma = s)))
#> [1] TRUE
```
