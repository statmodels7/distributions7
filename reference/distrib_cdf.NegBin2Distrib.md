# Negative Binomial Cumulative Distribution Function, NB2

Computes the negative binomial distribution function, the partial sum of
the mass \$\$F(q; \mu, \theta) = \sum\_{k=0}^{\lfloor q \rfloor} P(Y =
k; \mu, \theta),\$\$ by calling
[`stats::pnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at
`size = theta` and `mu = mu`, which evaluates it through the incomplete
beta function rather than by summing. The function is a step function,
constant between consecutive integers. Both tails are available exactly,
and `log.p = TRUE` returns a logarithm that stays finite where the
probability itself underflows.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

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

[`distrib_quantile.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.NegBin2Distrib.md)
for the generalized inverse,
[`distrib_pdf.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin2Distrib.md)
for the mass,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, which are an
exact finite sum for a discrete family, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
th <- list(mu = 4, theta = 2)

# The method is stats::pnbinom at size = theta, mu = mu.
all.equal(distrib_cdf(d, c(0, 2, 6), th), pnbinom(c(0, 2, 6), 2, mu = 4))
#> [1] TRUE

# A step function: it does not move between two consecutive integers.
distrib_cdf(d, c(2, 2.5, 2.999), th)
#> [1] 0.4074074 0.4074074 0.4074074

# It is the partial sum of the mass.
all.equal(distrib_cdf(d, 6, th), sum(distrib_pdf(d, 0:6, th)))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 6, th) + distrib_cdf(d, 6, th, lower.tail = FALSE)
#> [1] 1

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 2000, th, lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 2000, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -804.8314
```
