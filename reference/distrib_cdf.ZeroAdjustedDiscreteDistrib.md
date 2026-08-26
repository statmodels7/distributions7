# Zero-Adjusted Discrete Cumulative Distribution Function

Computes the hurdle distribution function \$\$F\_{ZA}(q) = \pi +
(1-\pi)\frac{F(q;\theta) - f(0;\theta)}{1 - f(0;\theta)} \quad (q \ge
0),\$\$ the parent's own distribution function shifted and rescaled by
the truncation. It is `0` below zero and is clamped to \\\[0, 1\]\\
against rounding.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- q:

  A numeric vector of quantiles. Values below zero give `0`.

- theta:

  A named list with the parent's parameters followed by `za`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, formed as \\1 - F\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, in \\\[0, 1\]\\.

## Notation

\\f\\ is the parent's mass function, \\F\\ its distribution function,
\\\pi\\ the probability of a zero, \\s\\ the parent's score and \\\ell\\
the log-mass of one observation.

## See also

[`distrib_pdf.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedDiscreteDistrib.md)
for the mass function,
[`distrib_quantile.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ZeroAdjustedDiscreteDistrib.md),
which inverts this, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)

distrib_cdf(d, c(0, 2, 5), theta)
#> [1] 0.4000000 0.6357806 0.9470111

# At zero it is the parameter itself, the whole mass there.
c(cdf_at_zero = distrib_cdf(d, 0, theta), parameter = 0.4)
#> cdf_at_zero   parameter 
#>         0.4         0.4 

# It agrees with the mass function summed, as it must on a lattice.
c(cdf = distrib_cdf(d, 5, theta), summed = sum(distrib_pdf(d, 0:5, theta)))
#>       cdf    summed 
#> 0.9470111 0.9470111 

# Both tails and the logarithm.
distrib_cdf(d, 2, theta, lower.tail = FALSE)
#> [1] 0.3642194
distrib_cdf(d, 2, theta, log.p = TRUE)
#> [1] -0.4529017
```
