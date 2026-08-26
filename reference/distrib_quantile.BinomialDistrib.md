# Binomial Quantile Function

Computes the generalized inverse of the binomial distribution function,
\$\$Q(p; \mu) = \min\\k \in \\0, \dots, n\\ : F(k; \mu) \ge p\\,\$\$ by
calling [`stats::qbinom()`](https://rdrr.io/r/stats/Binomial.html) at
`size = distrib@size`. The distribution function is a step function, so
the round trip through
[`distrib_cdf.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BinomialDistrib.md)
returns a probability **at least** `p`.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `p`. `mu` must lie in \\(0, 1)\\.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of integers between 0 and `size`, of length
`max(length(p), length(mu), length(distrib@size))`.

## See also

[`distrib_cdf.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BinomialDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)
th <- list(mu = 0.3)

# Integers, bounded by size.
distrib_quantile(d, c(0.025, 0.5, 0.975, 1), th)
#> [1]  0  3  6 10

# The round trip overshoots, the support being a lattice.
p <- c(0.025, 0.5, 0.975)
rbind(asked = p, reached = distrib_cdf(d, distrib_quantile(d, p, th), th))
#>               [,1]      [,2]      [,3]
#> asked   0.02500000 0.5000000 0.9750000
#> reached 0.02824752 0.6496107 0.9894079
```
