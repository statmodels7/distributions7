# Gumbel Random Number Generator

Draws `n` independent Gumbel variates as \\\mu - \sigma\log E\\ with
\\E\\ standard exponential, using
[`stats::rexp()`](https://rdrr.io/r/stats/Exponential.html). That is the
inverse-transform identity written without a logarithm of a uniform,
which keeps the left tail accurate: a uniform near zero has few
significant digits left after one logarithm and none after two.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one variate per parameter setting.
  `sigma` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws on the whole real line.

## See also

[`distrib_quantile.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GumbelDistrib.md)
for the inverse it rests on,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()

# The identity the method uses: mu - sigma log(E), E standard exponential.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 0, sigma = 1))
set.seed(2)
identical(a, 0 - 1 * log(rexp(3)))
#> [1] TRUE

# The sample moments recover the parameters through the fixed shape: the
# scale is sd sqrt(6)/pi and the location the mean less gamma times it.
set.seed(8)
z <- distrib_rng(d, 2e4, list(mu = 3, sigma = 2))
s <- sd(z) * sqrt(6) / pi
c(mu = mean(z) + digamma(1) * s, sigma = s)
#>       mu    sigma 
#> 3.005304 2.015630 
```
