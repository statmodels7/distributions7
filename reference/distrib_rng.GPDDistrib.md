# Generalized Pareto Random Generation

Draws by the inverse transform,
`distrib_quantile(distrib, runif(n), theta)`. The quantile function of
this family is elementary, so the transform is exact and costs one
uniform per draw; the base class's ratio-of-uniforms fallback is not
needed.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- n:

  A single positive integer, the number of draws.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of length `n`; a component of length 1 is recycled, so
  a parameter varying by observation draws one value per observation
  from its own member of the family.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws, non-negative and bounded above by
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md)
when the shape is negative.

## See also

[`distrib_quantile.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GPDDistrib.md),
which it inverts through,
[`distrib_pdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GPDDistrib.md)
for the density the draws follow, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
th <- list(sigma = 1.5, xi = 0.3)

set.seed(41)
x <- distrib_rng(d, 2e5, th)

# The mean exists here (xi < 1) and matches sigma/(1 - xi).
c(sample = mean(x), theory = mean(d, th))
#>   sample   theory 
#> 2.126313 2.142857 

# A negative shape draws inside a bounded support.
set.seed(42)
xn <- distrib_rng(d, 1e4, list(sigma = 2, xi = -0.4))
c(max_draw = max(xn), endpoint = distributions7:::gpd_endpoint(2, -0.4))
#> max_draw endpoint 
#> 4.845163 5.000000 

# Shape zero draws an exponential sample.
set.seed(43)
x0 <- distrib_rng(d, 1e5, list(sigma = 1.5, xi = 0))
c(sample_mean = mean(x0), theory = 1.5)
#> sample_mean      theory 
#>    1.491576    1.500000 
```
