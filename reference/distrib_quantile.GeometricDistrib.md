# Geometric Quantile Function

Computes the generalized inverse of the geometric distribution function,
\$\$Q(p; \mu) = \min\\k \in \\0, 1, 2, \dots\\ : F(k; \mu) \ge p\\,\$\$
by calling [`stats::qgeom()`](https://rdrr.io/r/stats/Geometric.html) at
`prob = 1/(1+mu)`. The distribution function is a step function, so the
round trip through
[`distrib_cdf.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GeometricDistrib.md)
returns a probability **at least** `p`.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

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

[`distrib_cdf.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GeometricDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
th <- list(mu = 3)

# Integers, and the median well below the mean, the mass being at zero.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1]  0  2 12

# The round trip overshoots, the support being a lattice.
p <- c(0.025, 0.5, 0.975)
rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
#>          [,1]     [,2]      [,3]
#> asked   0.025 0.500000 0.9750000
#> reached 0.250 0.578125 0.9762427
```
