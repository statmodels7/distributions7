# Mean of the Beta-Binomial Distribution

Closed form: \\E\[Y\] = n\mu\\, the same as a binomial's. Mixing the
success probability over a beta leaves the mean where it was, the beta's
own mean being \\\mu\\; the dispersion shows in the variance and above,
not here.

## Arguments

- x:

  A `BetaBinom1Distrib`, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with components `mu` (the success probability, strictly
  between 0 and 1) and `sigma` (the dispersion, positive), each a
  numeric vector of length 1 or `n`.

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

# n mu, and the dispersion does not move it.
mean(d, list(mu = 0.3, sigma = c(0.01, 0.5, 5)))
#> [1] 3 3 3
```
