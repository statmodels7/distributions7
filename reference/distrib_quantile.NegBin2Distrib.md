# Negative Binomial Quantile Function, NB2

Computes the generalized inverse of the distribution function, \$\$Q(p;
\mu, \theta) = \min\left\\y \in \mathbb{N}\_0 : F(y; \mu, \theta) \ge
p\right\\,\$\$ by calling
[`stats::qnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at
`size = theta` and `mu = mu`. The support is a lattice, so the round
trip through
[`distrib_cdf.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin2Distrib.md)
does **not** return `p`: it returns the mass up to the smallest integer
whose cumulative probability reaches `p`, which overshoots.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; `p = 1` gives `Inf`.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. Both must be strictly positive.

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
`max(length(p), length(mu), length(theta))`.

## See also

[`distrib_cdf.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin2Distrib.md),
which this inverts;
[`distrib_rng.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.NegBin2Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
th <- list(mu = 4, theta = 2)

# A nominal central 95 percent interval on the counts.
p <- c(0.025, 0.5, 0.975)
q <- distrib_quantile(d, p, th)
q
#> [1]  0  3 13

# On a lattice the round trip overshoots: the cumulative probability at the
# returned integer is at or above the one asked for, not equal to it.
rbind(asked = p, reached = distrib_cdf(d, q, th))
#>              [,1]      [,2]      [,3]
#> asked   0.0250000 0.5000000 0.9750000
#> reached 0.1111111 0.5390947 0.9805889

# The generalized inverse is a step function of p.
distrib_quantile(d, c(0.3, 0.4, 0.5), th)
#> [1] 2 2 3
```
