# NB1 Quantile Function

Computes the generalized inverse of the distribution function, \$\$Q(p;
\mu, \theta) = \min\left\\y \in \mathbb{N}\_0 : F(y; \mu, \theta) \ge
p\right\\,\$\$ by calling
[`stats::qnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at size
\\r = \mu/\theta\\ and success probability \\1/(1+\theta)\\. The support
is a lattice, so the round trip through
[`distrib_cdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin1Distrib.md)
does **not** return `p`: it returns the mass up to the smallest integer
whose cumulative probability reaches `p`, which overshoots.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

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

[`distrib_cdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin1Distrib.md),
which this inverts;
[`distrib_rng.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.NegBin1Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
th <- list(mu = 4, theta = 4)

# A nominal central 95 percent interval on the counts.
p <- c(0.025, 0.5, 0.975)
q <- distrib_quantile(d, p, th)
q
#> [1]  0  3 16

# On a lattice the round trip overshoots: the cumulative probability at the
# returned integer is at or above the one asked for, not equal to it.
rbind(asked = p, reached = distrib_cdf(d, q, th))
#>          [,1]   [,2]     [,3]
#> asked   0.025 0.5000 0.975000
#> reached 0.200 0.5904 0.977482

# The overdispersion shows in the width: the same mean under a Poisson
# reaches far less far.
c(nb1 = distrib_quantile(d, 0.975, th),
  poisson = distrib_quantile(poisson_distrib(), 0.975, list(mu = 4)))
#>     nb1 poisson 
#>      16       8 
```
