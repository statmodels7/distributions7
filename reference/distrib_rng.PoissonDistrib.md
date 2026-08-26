# Poisson Random Number Generator

Draws `n` independent Poisson counts by calling
[`stats::rpois()`](https://rdrr.io/r/stats/Poisson.html), so the draws
come from R's own generator and depend on `.Random.seed` in the usual
way. The cumulative-table fallback the discrete base class supplies is
bypassed.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of length `n`. A value of length 1 is recycled, so a
  vector of length `n` draws one count per mean. `mu` must be strictly
  positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An integer vector of `n` non-negative counts.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the mean back, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()

# Same generator as stats::rpois, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 5, list(mu = 3))
set.seed(2)
identical(a, rpois(5, lambda = 3))
#> [1] TRUE

# The sample mean and variance both estimate mu, this family being
# equidispersed.
set.seed(4)
z <- distrib_rng(d, 2e4, list(mu = 4.2))
c(mean = mean(z), var = var(z))
#>     mean      var 
#> 4.182950 4.165888 
```
