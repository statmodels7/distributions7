# Zero-Inflated Probability Mass Function

Computes the zero-inflated mass function \$\$P(Y = y) =
\zeta\\\mathbb{I}(y = 0) + (1-\zeta)\\ f(y; \theta),\$\$ with \\f\\ the
parent's mass function. Inflation can only ADD zeros: the mass at zero
is \\\zeta + (1-\zeta)f(0)\\, which exceeds \\f(0)\\ for every \\\zeta
\> 0\\, and no single observed zero can be attributed to one mechanism
or the other.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- y:

  A numeric vector of observations. A point off the parent's support
  returns whatever the parent returns there, scaled by \\1-\zeta\\.

- theta:

  A named list with the parent's parameters followed by `zi`, each a
  numeric vector of length 1 or of the length of `y`. `zi` must lie
  strictly inside \\(0, 1)\\.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned. The
  logarithm is taken of the mixture, not inside the parent, so it
  underflows where the mixture does. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the recycled length of `y` and `theta`.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score and
\\\ell\\ the log-mass of one observation.

## See also

[`distrib_cdf.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroInflatedDistrib.md)
for the distribution function,
[`distrib_gradient.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroInflatedDistrib.md)
for the score,
[`distrib_pdf.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedDiscreteDistrib.md)
for the wrapper that replaces the mass at zero, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)

distrib_pdf(d, 0:5, theta)
#> [1] 0.28734030 0.11202090 0.16803136 0.16803136 0.12602352 0.07561411

# Which is the mixture written out.
y <- 0:5
all.equal(distrib_pdf(d, y, theta), 0.25 * (y == 0) + 0.75 * dpois(y, 3))
#> [1] TRUE

# The mass at zero exceeds the parent's, inflation adding to it.
c(inflated = distrib_pdf(d, 0, theta), parent = dpois(0, 3))
#>   inflated     parent 
#> 0.28734030 0.04978707 

# And the whole mass function still sums to one.
sum(distrib_pdf(d, 0:200, theta))
#> [1] 1
```
