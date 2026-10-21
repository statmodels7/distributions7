# Excess Kurtosis of the Beta-Binomial Distribution

The method obtains the excess kurtosis from
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md)
as one quotient, so that it does not form a kurtosis near 3 and subtract
3. The success probability sets its sign. Near \\\mu = 1/2\\ it is
negative and falls towards \\-2\\ as the dispersion grows, because the
mass splits between the two endpoints and the law approaches a scaled
Bernoulli; near either end it is positive and large, because the long
tail towards the middle dominates.

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

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu \in (0,1)\\ is the success probability, \\\sigma \> 0\\ the
dispersion, \\n\\ the number of trials.

## See also

[`skewness.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom1Distrib.md),
from the same moments;
[`kurtosis.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BinomialDistrib.md),
the \\\sigma \to 0\\ limit;
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md);
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

## Examples

``` r
d <- betabinom1_distrib(size = 10)

# At a central success probability it falls towards -2 with the dispersion.
round(kurtosis(d, list(mu = 0.5, sigma = c(0.1, 1, 100))), 4)
#> [1] -0.5462 -1.5091 -1.9934

# At an extreme one it is positive and large.
round(kurtosis(d, list(mu = 0.05, sigma = c(0.01, 0.1, 1))), 4)
#> [1]  2.2225  7.2095 15.5435

# It falls onto the binomial's as the dispersion vanishes.
c(betabinomial = kurtosis(d, list(mu = 0.3, sigma = 1e-8)),
  binomial     = kurtosis(binomial_distrib(size = 10), list(mu = 0.3)))
#> betabinomial     binomial 
#>   -0.1238095   -0.1238095 
```
