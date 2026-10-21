# Excess Kurtosis of the Beta-Binomial Distribution in Two Shapes

The excess kurtosis is that of
[`kurtosis.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md)
at \\\mu = \alpha/(\alpha+\beta)\\ and \\\sigma = 1/(\alpha+\beta)\\,
computed by
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md).
Its sign follows the balance of the two shapes: negative where they are
comparable, positive where one dominates.

## Arguments

- x:

  A `BetaBinom2Distrib`, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `alpha` and `beta`, both positive; each
  is a numeric vector of length 1 or of the common length of the
  parameter settings.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$alpha), length(theta$beta))`.

## Notation

\\\alpha \> 0\\ and \\\beta \> 0\\ are the shapes of the mixing beta,
\\n\\ the number of trials.

## See also

[`skewness.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BetaBinom2Distrib.md),
from the same moments;
[`kurtosis.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BetaBinom1Distrib.md);
[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md);
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

## Examples

``` r
d <- betabinom2_distrib(size = 10)

# The excess kurtosis is negative where the shapes are comparable, positive
# where one dominates.
round(kurtosis(d, list(alpha = c(3, 0.5, 0.2), beta = 3)), 4)
#> [1] -0.7250  2.2219  9.3036

# The two parametrizations agree on the law.
all.equal(kurtosis(d, list(alpha = 2, beta = 3)),
          kurtosis(betabinom1_distrib(size = 10),
                   list(mu = 0.4, sigma = 1 / 5)))
#> [1] TRUE
```
