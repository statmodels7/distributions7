# Variance of the Beta-Binomial Distribution

Closed form: \$\$\operatorname{Var}(Y) = n\mu(1-\mu)\\\frac{1 +
n\sigma}{1 + \sigma}.\$\$ The first factor is the binomial variance and
the second is the overdispersion, at least 1 for every \\n \ge 1\\ and
growing towards \\n\\ as the dispersion rises. At \\\sigma \to 0\\ the
family is binomial; at large \\\sigma\\ the variance approaches
\\n^2\mu(1-\mu)\\, which is a Bernoulli scaled by the trial count and is
the largest a distribution on \\\\0, \ldots, n\\\\ with that mean can
have.

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

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu \in (0,1)\\ is the success probability, \\\sigma \> 0\\ the
dispersion and \\n\\ the number of trials.

## See also

[`mean.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom1Distrib.md);
[`variance.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.BinomialDistrib.md),
the \\\sigma \to 0\\ limit;
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md);
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

## Examples

``` r
d <- betabinom1_distrib(size = 10)

# The published form, written out.
all.equal(variance(d, list(mu = 0.3, sigma = 0.5)),
          10 * 0.3 * 0.7 * (1 + 10 * 0.5) / (1 + 0.5))
#> [1] TRUE

# It falls onto the binomial variance as the dispersion vanishes.
c(betabinomial = variance(d, list(mu = 0.3, sigma = 1e-8)),
  binomial     = variance(binomial_distrib(size = 10), list(mu = 0.3)))
#> betabinomial     binomial 
#>          2.1          2.1 
```
