# Skewness of the Poisson-Inverse Gaussian Distribution

Closed form from the cumulants, \\\gamma_1 = \kappa_3/\kappa_2^{3/2}\\
with \$\$\kappa_2 = \mu + \sigma\mu^2, \qquad \kappa_3 = \mu +
3\sigma\mu^2 + 3\sigma^2\mu^3.\$\$ This is where the family separates
from a negative binomial: matched on mean and variance, the
Poisson-inverse Gaussian is the more skewed of the two, its extra term
in \\\sigma^2\mu^3\\ having no counterpart there.

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

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma))`, positive throughout.

## Details

The cumulants of a mixed Poisson are the Poisson's plus what the mixing
distribution contributes, which is why they are the natural quantities
here and why a standardized moment follows from them without a sum over
the support.

## Notation

\\\mu \> 0\\ is the mean, \\\sigma \> 0\\ the dispersion and
\\\kappa_k\\ the \\k\\-th cumulant.

## See also

[`kurtosis.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Pig1Distrib.md),
from the fourth cumulant;
[`variance.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Pig1Distrib.md);
[`skewness.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin2Distrib.md),
the comparison;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

## Examples

``` r
d <- pig1_distrib()

# The cumulant ratio, written out.
all.equal(skewness(d, list(mu = 3, sigma = 1)),
          (3 + 3 * 9 + 3 * 27) / 12^1.5)
#> [1] TRUE

# Matched on mean and variance, it is more skewed than a negative binomial.
c(pig = skewness(d, list(mu = 3, sigma = 1)),
  negbin = skewness(negbin2_distrib(), list(mu = 3, theta = 1)))
#>      pig   negbin 
#> 2.670245 2.020726 
```
