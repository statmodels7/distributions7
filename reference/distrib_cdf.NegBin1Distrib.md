# NB1 Cumulative Distribution Function

Computes the negative binomial distribution function, the partial sum of
the mass, by calling
[`stats::pnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at size
\\r = \mu/\theta\\ and success probability \\1/(1+\theta)\\. That
function evaluates it through the incomplete beta function, so nothing
is summed. The result is a step function, constant between consecutive
integers. Both tails are available exactly, and `log.p = TRUE` returns a
logarithm that stays finite where the probability itself underflows.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- q:

  A numeric vector of quantiles. A non-integer value is floored, and a
  value below zero gives a lower-tail probability of 0.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. Both must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu), length(theta))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBin1Distrib.md)
for the generalized inverse,
[`distrib_pdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md)
for the mass,
[`nb1_size()`](https://statmodels7.github.io/distributions7/reference/nb1_size.md)
and
[`nb1_prob()`](https://statmodels7.github.io/distributions7/reference/nb1_prob.md)
for the pairing this uses, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
th <- list(mu = 4, theta = 4)

# The method is stats::pnbinom at size mu/theta and prob 1/(1 + theta).
all.equal(distrib_cdf(d, c(0, 2, 6), th),
          pnbinom(c(0, 2, 6), size = 1, prob = 1 / 5))
#> [1] TRUE

# A step function: it does not move between two consecutive integers.
distrib_cdf(d, c(2, 2.5, 2.999), th)
#> [1] 0.488 0.488 0.488

# It is the partial sum of the mass.
all.equal(distrib_cdf(d, 6, th), sum(distrib_pdf(d, 0:6, th)))
#> [1] TRUE

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 2000, th, lower.tail = FALSE)
#> [1] 1.210776e-194
distrib_cdf(d, 2000, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -446.5102
```
