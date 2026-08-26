# Excess Kurtosis of the Poisson-Inverse Gaussian Distribution

Closed form from the cumulants, \\\gamma_2 = \kappa_4/\kappa_2^{2}\\
with \$\$\kappa_2 = \mu + \sigma\mu^2, \qquad \kappa_4 = \mu +
7\sigma\mu^2 + 18\sigma^2\mu^3 + 15\sigma^3\mu^4.\$\$ The ratio of the
fourth cumulant to the squared second **is** the excess kurtosis, so no
3 is subtracted here: the subtraction is what turns a fourth *moment*
into an excess, and a fourth cumulant already carries it.

## Arguments

- x:

  A `Pig1Distrib`, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `sigma` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$sigma))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean, \\\sigma \> 0\\ the dispersion and
\\\kappa_k\\ the \\k\\-th cumulant.

## See also

[`skewness.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Pig1Distrib.md),
from the third cumulant;
[`kurtosis.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin2Distrib.md),
the comparison;
[`kurtosis.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PoissonDistrib.md),
the limit;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

## Examples

``` r
d <- pig1_distrib()

# The cumulant ratio, written out.
all.equal(kurtosis(d, list(mu = 3, sigma = 1)),
          (3 + 7 * 9 + 18 * 27 + 15 * 81) / 144)
#> [1] TRUE

# It falls onto the Poisson's 1 / mu as the dispersion vanishes.
c(pig = kurtosis(d, list(mu = 3, sigma = 1e-9)),
  poisson = kurtosis(poisson_distrib(), list(mu = 3)))
#>       pig   poisson 
#> 0.3333333 0.3333333 
```
