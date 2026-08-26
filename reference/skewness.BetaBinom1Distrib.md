# Skewness of the Beta-Binomial Distribution

Assembled from the falling factorial moments, \\\gamma_1 =
c_3/c_2^{3/2}\\ with \\c_k\\ the central moments
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md)
returns. It is zero at \\\mu = 1/2\\ whatever the dispersion, the mixing
beta being symmetric there, and it takes the sign of \\1 - 2\mu\\
elsewhere.

## Arguments

- x:

  A `BetaBinom1Distrib`, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `mu` (strictly between 0 and 1) and
  `sigma` (positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma))`.

## Details

The falling factorial moments are the quantities a beta-binomial has in
closed form: each is a product of \\k\\ ratios and needs no sum over the
support. The raw and then the central moments follow from them by two
written-out linear maps.

## Notation

\\\mu \in (0,1)\\ is the success probability, \\\sigma \> 0\\ the
dispersion, \\n\\ the number of trials and \\c_k\\ the \\k\\-th central
moment.

## See also

[`kurtosis.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md),
from the same moments;
[`skewness.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BinomialDistrib.md),
the \\\sigma \to 0\\ limit;
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md);
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

## Examples

``` r
d <- betabinom1_distrib(size = 10)

# Zero at a success probability of one half, at any dispersion.
round(skewness(d, list(mu = 0.5, sigma = c(0.1, 1, 10))), 12)
#> [1] 0 0 0

# It falls onto the binomial's as the dispersion vanishes.
c(betabinomial = skewness(d, list(mu = 0.3, sigma = 1e-8)),
  binomial     = skewness(binomial_distrib(size = 10), list(mu = 0.3)))
#> betabinomial     binomial 
#>    0.2760263    0.2760262 
```
