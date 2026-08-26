# Variance of the Beta-Binomial Distribution in Two Shapes

The second central moment, assembled from the falling factorial moments
by
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md)
at \\\mu = \alpha/(\alpha+\beta)\\ and \\\sigma = 1/(\alpha+\beta)\\. It
exceeds the binomial variance at the same mean, and the excess grows as
the two shapes shrink. That excess is the family's reason for being.

## Arguments

- x:

  A `BetaBinom2Distrib`, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `alpha` and `beta`, both positive, each a
  numeric vector of length 1 or `n`. Large shapes concentrate the mixing
  beta and take the family towards a binomial.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the shapes of the mixing beta and
\\n\\ the number of trials.

## See also

[`mean.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom2Distrib.md),
[`skewness.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom2Distrib.md);
[`variance.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom1Distrib.md);
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md);
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

## Examples

``` r
d <- betabinom2_distrib(size = 10)

# The two parametrizations agree on the law.
all.equal(variance(d, list(alpha = 2, beta = 3)),
          variance(betabinom1_distrib(size = 10),
                   list(mu = 0.4, sigma = 1 / 5)))
#> [1] TRUE

# It exceeds the binomial variance at the same mean.
c(betabinomial = variance(d, list(alpha = 2, beta = 3)),
  binomial     = variance(binomial_distrib(size = 10), list(mu = 0.4)))
#> betabinomial     binomial 
#>          6.0          2.4 
```
