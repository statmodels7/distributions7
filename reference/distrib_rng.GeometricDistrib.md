# Geometric Random Number Generator

Draws `n` independent geometric counts by calling
[`stats::rgeom()`](https://rdrr.io/r/stats/Geometric.html) at
`prob = 1/(1+mu)`, so the draws come from R's own generator and depend
on `.Random.seed` in the usual way.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

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
d <- geometric_distrib()

# Same generator as stats::rgeom, so the same seed gives the same draws.
set.seed(2)
a <- distrib_rng(d, 5, list(mu = 3))
set.seed(2)
identical(a, rgeom(5, prob = 1 / 4))
#> [1] TRUE

# The sample mean estimates mu and the variance mu(1+mu).
set.seed(5)
z <- distrib_rng(d, 2e4, list(mu = 2.5))
c(mean = mean(z), var = var(z), mu_1_plus_mu = 2.5 * 3.5)
#>         mean          var mu_1_plus_mu 
#>      2.48200      8.68201      8.75000 
```
