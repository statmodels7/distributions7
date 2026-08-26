# Folded Random Generation

Draws `n` values from the parent and takes the absolute value. That is
the DEFINITION of the folded variable, not an approximation of it, so
the draws are exact whatever route the parent's own generator takes.
They consume whatever the parent consumes from R's stream, and are
reproducible under
[`base::set.seed()`](https://rdrr.io/r/base/Random.html).

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- n:

  The number of draws, a single non-negative whole number.

- theta:

  A named list of the parent's parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `n`, every value non-negative.

## See also

[`distrib_pdf.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.FoldedDistrib.md)
for the density these are drawn from,
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
for the family, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)

set.seed(1)
distrib_rng(d, 5, theta)
#> [1] 0.2517446 0.7203720 0.5027543 2.4143370 0.8954093

# It is the parent's draw with the sign removed, from the same seed.
set.seed(1)
abs(distrib_rng(gaussian1_distrib(), 5, theta))
#> [1] 0.2517446 0.7203720 0.5027543 2.4143370 0.8954093

# And a large sample reproduces the folded density.
set.seed(2)
big <- distrib_rng(d, 20000, theta)
c(sampled = mean(big < 1),
  exact = distrib_cdf(d, 1, theta))
#>   sampled     exact 
#> 0.5517000 0.5558891 
```
