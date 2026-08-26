# Zero-Inflated Cumulative Distribution Function

Computes \$\$F\_{ZI}(q) = (1-\zeta) F(q; \theta) + \zeta\\\mathbb{I}(q
\ge 0)\$\$ from the parent's own distribution function, exactly and with
no summation. The result is clamped to \\\[0, 1\]\\ against rounding.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- q:

  A numeric vector of quantiles. Values below zero give `0`.

- theta:

  A named list with the parent's parameters followed by `zi`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, formed as \\1 - F\\ and
  so subject to that subtraction's cancellation far into the upper tail.

- log.p:

  Logical of length 1. When `TRUE` the logarithm is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, in \\\[0, 1\]\\.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score and
\\\ell\\ the log-mass of one observation.

## See also

[`distrib_pdf.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroInflatedDistrib.md)
for the mass function,
[`distrib_quantile.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroInflatedDistrib.md),
which inverts this, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)

q <- c(0, 2, 5)
distrib_cdf(d, q, theta)
#> [1] 0.2873403 0.5673926 0.9370615

# Which is the parent's, shrunk and shifted.
all.equal(distrib_cdf(d, q, theta), 0.75 * ppois(q, 3) + 0.25)
#> [1] TRUE

# It agrees with the mass function summed, as it must on a lattice.
c(cdf = distrib_cdf(d, 5, theta), summed = sum(distrib_pdf(d, 0:5, theta)))
#>       cdf    summed 
#> 0.9370615 0.9370615 

# Both tails and the logarithm.
distrib_cdf(d, 2, theta, lower.tail = FALSE)
#> [1] 0.4326074
distrib_cdf(d, 2, theta, log.p = TRUE)
#> [1] -0.5667039
```
