# Skewness of the Beta-Binomial Distribution

The skewness has the closed form \$\$\gamma_1 = (1 - 2\mu)\\\frac{1 +
2n\sigma}{1 + 2\sigma} \sqrt{\frac{1 + \sigma}{n\mu(1-\mu)(1 +
n\sigma)}},\$\$ computed by
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md).
It is zero at \\\mu = 1/2\\ whatever the dispersion, because the mixing
beta is symmetric there, and it takes the sign of \\1 - 2\mu\\
elsewhere. The closed form forms no difference of raw moments, so it
keeps its accuracy at a small dispersion and a large \\n\\.

## Arguments

- x:

  A `BetaBinom1Distrib`, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `mu` (strictly between 0 and 1) and
  `sigma` (positive); each is a numeric vector of length 1 or of the
  common length of the parameter settings.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu \in (0,1)\\ is the success probability, \\\sigma \> 0\\ the
dispersion, \\n\\ the number of trials.

## See also

[`kurtosis.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md),
from the same moments;
[`skewness.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BinomialDistrib.md),
the \\\sigma \to 0\\ limit;
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md);
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

## Examples

``` r
d <- betabinom1_distrib(size = 10)

# The skewness is zero at a success probability of one half, at any
# dispersion.
round(skewness(d, list(mu = 0.5, sigma = c(0.1, 1, 10))), 12)
#> [1] 0 0 0

# It falls onto the binomial's as the dispersion vanishes.
c(betabinomial = skewness(d, list(mu = 0.3, sigma = 1e-8)),
  binomial     = skewness(binomial_distrib(size = 10), list(mu = 0.3)))
#> betabinomial     binomial 
#>    0.2760263    0.2760262 
```
