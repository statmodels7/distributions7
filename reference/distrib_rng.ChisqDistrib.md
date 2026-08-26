# Chi-Squared Random Number Generator

Draws `n` independent chi-squared variates by calling
[`stats::rchisq()`](https://rdrr.io/r/stats/Chisquare.html) at
`df = mu`, so the draws come from R's own generator and depend on
`.Random.seed` in the usual way. The ratio-of-uniforms fallback the base
class supplies is bypassed.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of length `n`. A component of length 1 is recycled, so a vector of
  length `n` draws one variate per parameter setting. It must be
  strictly positive and need not be a whole number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` strictly positive draws.

## See also

[`distrib_quantile.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ChisqDistrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameter back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()

# Same generator as stats::rchisq at df = mu.
set.seed(2)
a <- distrib_rng(d, 3, list(mu = 4))
set.seed(2)
identical(a, rchisq(3, df = 4))
#> [1] TRUE

# The sample mean recovers the degrees of freedom, and the variance is
# about twice it, the two being tied.
set.seed(7)
z <- distrib_rng(d, 2e4, list(mu = 4))
c(mean = mean(z), var = var(z))
#>     mean      var 
#> 3.997527 7.986421 
```
