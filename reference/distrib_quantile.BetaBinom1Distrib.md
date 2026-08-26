# Beta-Binomial Quantile Function

Computes the generalized inverse \\Q(p) = \min\\y : F(y) \ge p\\\\ by
walking the exact cumulative sum over the support. The result is an
integer count, the distribution function being a step function, so
`Q(F(y))` returns `y` while `F(Q(p))` is at least `p` and generally
exceeds it. A tolerance of `1e-12` is allowed on the comparison so that
a `p` obtained from
[`distrib_cdf.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BetaBinom1Distrib.md)
maps back to the count it came from.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. `NA` is returned for `NA`; a value at
  or below 0 gives 0 and one at or above 1 gives `n`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1. `mu` must lie in \\(0, 1)\\ and `sigma` be strictly
  positive. One cumulative table is built for the whole call, so a
  parameter varying by observation is not supported here.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is read as a logarithm. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of counts in \\\\0, \dots, n\\\\, of length
`length(p)`.

## See also

[`distrib_cdf.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BetaBinom1Distrib.md)
for the function inverted here,
[`distrib_rng.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BetaBinom1Distrib.md)
for draws, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)

# The deciles, which are counts.
distrib_quantile(d, c(0.1, 0.5, 0.9), th)
#> [1] 0 2 8

# The round trip from a count through F and back is the identity.
all.equal(distrib_quantile(d, distrib_cdf(d, 0:10, th), th), as.numeric(0:10))
#> [1] TRUE

# The other direction only reaches at least p, F being a step function.
rbind(p = c(0.2, 0.45, 0.8),
      F_at_Q = distrib_cdf(d, distrib_quantile(d, c(0.2, 0.45, 0.8), th), th))
#>             [,1]      [,2]      [,3]
#> p      0.2000000 0.4500000 0.8000000
#> F_at_Q 0.2644607 0.5338989 0.8445202
```
