# Zero-Adjusted Continuous Probability Density Function

Computes the mixed density \$\$f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi)
f_W(y;\theta) \quad (y \ne 0).\$\$ The parent needs no truncation,
because a continuous parent has \\P(W = 0) = 0\\: the atom can be placed
without removing any mass, and the density is simply scaled by
\\1-\pi\\.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations. Exactly zero gives \\\pi\\; any
  other value gives \\(1-\pi)\\ times the parent's density there,
  including a point outside the parent's support, where the parent's
  density is zero.

- theta:

  A named list with the parent's parameters followed by `za`; each
  element is a numeric vector of length 1 or of the length of `y`.

- log:

  Logical of length 1. When `TRUE` the logarithm is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the recycled length of `y` and `theta`.

## Details

The value at zero is a PROBABILITY and the values elsewhere are
DENSITIES, so the returned vector mixes two kinds of number. This is the
form of a mixed distribution, and the object therefore declares
[`distrib_atoms.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md):
without that declaration a consumer would integrate the returned
function and find \\1-\pi\\.

## Notation

\\f_W\\ is the parent's density, \\\pi\\ the probability of the atom at
zero, \\f_Y\\ the mixed density and \\\ell\\ the log-density of one
observation.

## See also

[`distrib_atoms.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md)
for the declaration that makes this a mixed distribution,
[`distrib_cdf.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroAdjustedContinuousDistrib.md)
for the distribution function, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

# The value at zero is a probability, and the values elsewhere are densities.
c(at_zero = distrib_pdf(d, 0, theta),
  at_two = distrib_pdf(d, 2, theta),
  scaled_parent = 0.7 * dnorm(2, 1, 2))
#>       at_zero        at_two scaled_parent 
#>     0.3000000     0.1232229     0.1232229 

# The density part integrates to 1 - pi, and the atom carries the rest.
integrate(function(z) ifelse(z == 0, 0, distrib_pdf(d, z, theta)),
          -Inf, Inf)$value
#> [1] 0.7

# For this reason the object declares its atom.
distrib_atoms(d, theta)
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 
```
