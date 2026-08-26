# Skew Normal Random Generation in the Centered Parametrization

Draws from the skew normal at the direct parameters the centered triple
implies, through
[`distrib_rng.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal1Distrib.md)'s
stochastic representation. The draws are exact and cost two `rnorm`
calls, whatever the skewness is.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, each a
  numeric vector of length 1 or of length `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws, whose first three sample moments estimate
`mu`, `sigma^2` and `gamma1`.

## See also

[`distrib_rng.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.SkewNormal1Distrib.md)
for the representation,
[`mean.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal2Distrib.md)
for the moments the draws reproduce, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
th <- list(mu = 3, sigma = 2, gamma1 = 0.6)

set.seed(2)
x <- distrib_rng(d, 2e5, th)

# The three parameters are the three sample moments, which is the whole
# point of this parametrization.
rbind(sample = c(mean(x), sd(x), mean((x - mean(x))^3) / sd(x)^3),
      parameter = c(th$mu, th$sigma, th$gamma1))
#>               [,1]     [,2]      [,3]
#> sample    2.997335 2.001777 0.5901477
#> parameter 3.000000 2.000000 0.6000000
```
