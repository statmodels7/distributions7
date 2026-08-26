# Negative Binomial Random Number Generator, NB2

Draws `n` independent negative binomial counts by calling
[`stats::rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at
`size = theta` and `mu = mu`, so the draws come from R's own generator
and depend on `.Random.seed` in the usual way. The cumulative-table
fallback the base class supplies for a discrete family is bypassed.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

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

[`distrib_quantile.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBin2Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()

# Same generator as stats::rnbinom at size = theta, mu = mu.
set.seed(2)
a <- distrib_rng(d, 5, list(mu = 4, theta = 2))
set.seed(2)
identical(a, rnbinom(5, size = 2, mu = 4))
#> [1] TRUE

# The sample moments recover the parameters: the mean directly, and the
# dispersion as mu^2/(var - mu).
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 4, theta = 2))
c(mu = mean(z), theta = mean(z)^2 / (var(z) - mean(z)))
#>       mu    theta 
#> 4.006700 2.041968 

# Overdispersion is the point: the variance is three times the mean here.
c(mean = mean(z), var = var(z))
#>     mean      var 
#>  4.00670 11.86855 
```
