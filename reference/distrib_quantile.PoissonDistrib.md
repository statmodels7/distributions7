# Poisson Quantile Function

Computes the generalized inverse of the Poisson distribution function,
\$\$Q(p; \mu) = \min\\k \in \\0, 1, 2, \dots\\ : F(k; \mu) \ge p\\,\$\$
by calling [`stats::qpois()`](https://rdrr.io/r/stats/Poisson.html). The
distribution function is a step function, so its inverse is a step
function too and the round trip through
[`distrib_cdf.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PoissonDistrib.md)
returns a probability **at least** `p`, not `p` itself.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. `p = 1` gives `Inf`.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `p`. `mu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of non-negative integers, of length
`max(length(p), length(mu))`.

## See also

[`distrib_cdf.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PoissonDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
th <- list(mu = 3)

# Integers, being the smallest count whose cumulative mass reaches p.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1] 0 3 7

# The round trip overshoots, the support being a lattice.
p <- c(0.025, 0.5, 0.975)
rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
#>               [,1]      [,2]      [,3]
#> asked   0.02500000 0.5000000 0.9750000
#> reached 0.04978707 0.6472319 0.9880955

# It is the smallest k with F(k) >= p, which the definition checks.
k <- distrib_quantile(d, 0.9, th)
c(k = k, F_k = distrib_cdf(d, k, th), F_k_minus_1 = distrib_cdf(d, k - 1, th))
#>           k         F_k F_k_minus_1 
#>   5.0000000   0.9160821   0.8152632 
```
