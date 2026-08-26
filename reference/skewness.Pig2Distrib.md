# Skewness of the Orthogonal Poisson-Inverse Gaussian Distribution

The cumulant ratio \\\gamma_1 = \kappa_3/\kappa_2^{3/2}\\ of
[`skewness.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig1Distrib.md),
read at the dispersion \\\sigma = \\
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)`(mu, alpha)`.
A standardized moment is a property of the law, so it agrees with the
other parametrization once the mapping is applied.

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

A numeric vector, of length
`max(length(theta$mu), length(theta$alpha))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean, \\\alpha \> 0\\ the orthogonal dispersion and
\\\kappa_k\\ the \\k\\-th cumulant.

## See also

[`kurtosis.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Pig2Distrib.md),
from the fourth cumulant;
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md);
[`skewness.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig1Distrib.md);
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

## Examples

``` r
d <- pig2_distrib()

# The two parametrizations agree once the mapping is applied.
s <- distributions7:::pig2_sigma(3, 1)
all.equal(skewness(d, list(mu = 3, alpha = 1)),
          skewness(pig1_distrib(), list(mu = 3, sigma = s)))
#> [1] TRUE
```
