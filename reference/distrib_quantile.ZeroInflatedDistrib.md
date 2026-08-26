# Zero-Inflated Quantile Function

Inverts the mixture distribution function. The quantile is `0` for \\p
\le \zeta + (1-\zeta)F(0; \theta)\\, which is the whole of the inflated
mass at zero, and otherwise the parent's quantile at the rescaled
probability \\(p - \zeta)/(1 - \zeta)\\.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- p:

  A numeric vector of probabilities, clamped to \\\[0, 1\]\\ after the
  `log.p` and `lower.tail` transformations are applied.

- theta:

  A named list with the parent's parameters followed by `zi`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm.
  Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of the recycled length of `p` and
`theta`.

## Details

The result is a LATTICE quantile and overshoots, as any discrete
quantile does: `distrib_cdf(d, distrib_quantile(d, p, theta), theta)`
returns the smallest attainable probability at or above `p`, not `p`
itself. The inflated atom makes the first step larger than the parent's,
so a probability below \\L_0\\ maps to zero however small it is.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score and
\\\ell\\ the log-mass of one observation.

## See also

[`distrib_cdf.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ZeroInflatedDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)

distrib_quantile(d, c(0.1, 0.3, 0.9), theta)
#> [1] 0 1 5

# Everything below the inflated mass at zero maps to zero.
c(mass_at_zero = distrib_pdf(d, 0, theta),
  q_below = distrib_quantile(d, 0.2, theta))
#> mass_at_zero      q_below 
#>    0.2873403    0.0000000 

# The round trip overshoots, as on any lattice: asked for 0.3, the
# attainable probability at the returned point is higher.
p <- c(0.1, 0.3, 0.9)
rbind(asked = p,
      reached = distrib_cdf(d, distrib_quantile(d, p, theta), theta))
#>              [,1]      [,2]      [,3]
#> asked   0.1000000 0.3000000 0.9000000
#> reached 0.2873403 0.3993612 0.9370615
```
