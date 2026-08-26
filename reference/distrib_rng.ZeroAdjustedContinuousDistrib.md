# Zero-Adjusted Continuous Random Number Generator

Draws zeros with probability \\\pi\\ and otherwise samples from the
parent unchanged. No truncation is involved, unlike the discrete branch:
a continuous parent produces an exact zero with probability zero, so a
draw from it never collides with the atom.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- n:

  The number of draws, a single non-negative whole number.

- theta:

  A named list with the parent's parameters followed by `za`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `n`, in which the value `0` appears with
probability \\\pi\\ and is exact.

## Notation

\\f_W\\ is the parent's density, \\\pi\\ the probability of the atom at
zero, \\f_Y\\ the mixed density and \\\ell\\ the log-density of one
observation.

## See also

[`distrib_pdf.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedContinuousDistrib.md)
for the law these are drawn from,
[`distrib_rng.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.ZeroAdjustedDiscreteDistrib.md),
which does truncate, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

set.seed(3)
distrib_rng(d, 8, theta)
#> [1]  0.000000  1.391566  1.060248  1.170835  3.233220 -1.437715  0.000000
#> [8]  0.000000

# A large sample reproduces the atom, and the non-zero draws reproduce the
# parent.
set.seed(3)
big <- distrib_rng(d, 20000, theta)
c(sampled_atom = mean(big == 0), parameter = 0.3)
#> sampled_atom    parameter 
#>      0.30195      0.30000 
round(c(mean = mean(big[big != 0]), sd = sd(big[big != 0])), 2)
#> mean   sd 
#> 0.97 2.03 
```
