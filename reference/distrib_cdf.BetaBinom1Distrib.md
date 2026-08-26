# Beta-Binomial Cumulative Distribution Function

Computes \\F(q) = P(Y \le q) = \sum\_{k \le \lfloor q \rfloor} P(Y =
k)\\ as the cumulative sum of the mass over the whole support \\\\0,
\dots, n\\\\. The support being finite, that sum is exact rather than an
approximation, and it costs \\n+1\\ evaluations of the mass whatever the
length of `q`. Below the support the probability is 0 and at or above
\\n\\ it is 1.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- q:

  A numeric vector of quantiles. A non-integer value is floored, so
  `F(2.7)` is `F(2)`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1. `mu` must lie in \\(0, 1)\\ and `sigma` be strictly
  positive. A parameter varying by observation is not supported here:
  one cumulative table is built for the whole call.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, formed as \\1 - F\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`length(q)`. With `log.p = TRUE` the values are logarithms and are
non-positive.

## See also

[`distrib_quantile.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BetaBinom1Distrib.md)
for the generalized inverse,
[`distrib_pdf.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom1Distrib.md)
for the mass, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)

# The cumulative sum of the mass, exactly.
all.equal(distrib_cdf(d, 0:10, th), cumsum(distrib_pdf(d, 0:10, th)))
#> [1] TRUE

# Below the support it is zero and at the top of it one.
distrib_cdf(d, c(-1, 0, 3, 10, 11), th)
#> [1] 0.0000000 0.2644607 0.6303587 1.0000000 1.0000000

# A non-integer quantile is floored.
c(distrib_cdf(d, 2.7, th), distrib_cdf(d, 2, th))
#> [1] 0.5338989 0.5338989

# The two tails sum to one.
distrib_cdf(d, 4, th) + distrib_cdf(d, 4, th, lower.tail = FALSE)
#> [1] 1
```
