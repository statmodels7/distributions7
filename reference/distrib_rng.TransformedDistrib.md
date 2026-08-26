# Transformed Random Number Generator

Draws `n` values from the parent and applies \\g\\ to them. That is the
DEFINITION of \\Y = g(X)\\, not an approximation of it, so the draws are
exact whatever route the parent's own generator takes and consume
exactly what the parent consumes from R's stream.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- n:

  The number of draws, a single non-negative whole number.

- theta:

  A named list of the parent's parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `n`, on the transformed scale.

## See also

[`distrib_pdf.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TransformedDistrib.md)
for the density these are drawn from, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)

set.seed(1)
distrib_rng(d, 5, theta)
#> [1] 0.9988376 1.9096398 0.8449288 5.9075113 2.1460012

# It is the parent's draw mapped forward, from the same seed.
set.seed(1)
exp(distrib_rng(gaussian1_distrib(), 5, theta))
#> [1] 0.9988376 1.9096398 0.8449288 5.9075113 2.1460012

# And a large sample reproduces the transformed distribution function.
set.seed(2)
big <- distrib_rng(d, 20000, theta)
c(sampled = mean(big < 2), exact = distrib_cdf(d, 2, theta))
#>   sampled     exact 
#> 0.5911500 0.5953906 
```
