# Mean of the Beta-Binomial Distribution

The mean has the closed form \\E\[Y\] = n\mu\\, the same as a
binomial's. The mixing beta has mean \\\mu\\, so mixing the success
probability over it does not change the mean; the dispersion enters the
variance and the higher moments.

## Arguments

- x:

  A `BetaBinom1Distrib`, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `mu` (the success probability, strictly
  between 0 and 1) and `sigma` (the dispersion, positive); each is a
  numeric vector of length 1 or of the common length of the parameter
  settings.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu \in (0,1)\\ is the success probability, \\\sigma \> 0\\ the
dispersion and \\n\\ the number of trials, held on the object.

## See also

[`variance.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom1Distrib.md),
where the dispersion does enter;
[`mean.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.BinomialDistrib.md),
which this equals;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

## Examples

``` r
d <- betabinom1_distrib(size = 10)

# The mean is n mu, and the dispersion does not move it.
mean(d, list(mu = 0.3, sigma = c(0.01, 0.5, 5)))
#> [1] 3 3 3
```
