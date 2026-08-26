# NB1 Random Number Generator

Draws `n` independent counts by calling
[`stats::rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at size
\\r = \mu/\theta\\ and success probability \\1/(1+\theta)\\, so the
draws come from R's own generator and depend on `.Random.seed` in the
usual way. The cumulative-table fallback the base class supplies for a
discrete family is bypassed.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled, so
  a vector of length `n` draws one count per parameter setting. Both
  must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` non-negative integers.

## See also

[`distrib_quantile.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBin1Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()

# Same generator as stats::rnbinom at size mu/theta and prob 1/(1 + theta).
set.seed(2)
a <- distrib_rng(d, 5, list(mu = 4, theta = 4))
set.seed(2)
identical(a, rnbinom(5, size = 1, prob = 1 / 5))
#> [1] TRUE

# The sample moments recover the parameters: the mean directly, and the
# dispersion as var/mean - 1.
set.seed(5)
z <- distrib_rng(d, 2e4, list(mu = 4, theta = 4))
c(mu = mean(z), theta = var(z) / mean(z) - 1)
#>       mu    theta 
#> 4.002300 4.126058 
```
