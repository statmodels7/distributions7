# Pseudo-Huber Random Number Generator

Draws `n` independent variates by inverse transform: uniform variates
from [`stats::runif()`](https://rdrr.io/r/stats/Uniform.html) passed
through
[`distrib_quantile.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md).
Each draw therefore costs a root-find over a quadrature, which makes
this the slowest generator in the package; a sample of a few thousand is
comfortable, a sample of a million is not. The draws depend on
`.Random.seed` in the usual way.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of length `n`. A component of length 1 is
  recycled, so a vector of length `n` draws one variate per parameter
  setting. `sigma` and `nu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## See also

[`distrib_quantile.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md)
for the inversion,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The draws are the quantile function at uniform variates, which is the
# whole mechanism and the whole cost.
set.seed(6)
a <- distrib_rng(d, 5, th)
set.seed(6)
identical(a, distrib_quantile(d, runif(5), th))
#> [1] TRUE

# A sample of a few hundred is comfortable, and its moments sit where the
# sampling error of that size puts them.
set.seed(6)
z <- distrib_rng(d, 300, th)
rbind(sample = c(mean(z), var(z)),
      theoretical = c(mean(d, th), variance(d, th)))
#>                  [,1]     [,2]
#> sample      0.6215911 5.079847
#> theoretical 0.4000000 4.429997
```
