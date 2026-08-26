# Skewness of the Beta-Binomial Distribution in Two Shapes

The standardized third central moment, \\\gamma_1 = c_3/c_2^{3/2}\\ with
the central moments \\c_k\\ from
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md)
at \\\mu = \alpha/(\alpha+\beta)\\ and \\\sigma = 1/(\alpha+\beta)\\. It
is zero when the two shapes are equal, the mixing beta being symmetric
there.

## Arguments

- x:

  A `BetaBinom2Distrib`, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `alpha` and `beta`, both positive, each a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the shapes of the mixing beta,
\\n\\ the number of trials and \\c_k\\ the \\k\\-th central moment.

## See also

[`kurtosis.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom2Distrib.md),
from the same moments;
[`skewness.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom1Distrib.md);
[`betabinom_central()`](https://statmodels7.github.io/distributions7/reference/betabinom_central.md);
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

## Examples

``` r
d <- betabinom2_distrib(size = 10)

# Zero when the two shapes are equal.
round(skewness(d, list(alpha = 3, beta = 3)), 12)
#> [1] 0

# The two parametrizations agree on the law.
all.equal(skewness(d, list(alpha = 2, beta = 3)),
          skewness(betabinom1_distrib(size = 10),
                   list(mu = 0.4, sigma = 1 / 5)))
#> [1] TRUE
```
