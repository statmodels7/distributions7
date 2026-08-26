# Weibull Cumulative Distribution Function

Computes the Weibull distribution function \$\$F(q; \mu, \sigma) = 1 -
\exp\left\\-(q/\mu)^{\sigma}\right\\\$\$ by calling
[`stats::pweibull()`](https://rdrr.io/r/stats/Weibull.html). The
survival function is elementary, \\1 - F(q) =
\exp\\-(q/\mu)^{\sigma}\\\\, so `lower.tail = FALSE` with `log.p = TRUE`
returns \\-(q/\mu)^{\sigma}\\ exactly and never forms the difference.
That combination is what a right-censored observation far out in the
tail contributes to a log-likelihood.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- q:

  A numeric vector of quantiles. Below the support the probability is 0
  in the lower tail and 1 in the upper.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
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
`max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Weibull1Distrib.md)
for the inverse,
[`distrib_pdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull1Distrib.md)
for the density,
[`distrib_grad_cdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Weibull1Distrib.md)
for the derivatives of this function in the parameters, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
th <- list(mu = 2, sigma = 1.5)

# The method is stats::pweibull at this parametrization.
all.equal(distrib_cdf(d, c(0.5, 1.2, 3.0), th),
          pweibull(c(0.5, 1.2, 3.0), shape = 1.5, scale = 2))
#> [1] TRUE

# The log survival function is the exponent itself, with no cancellation.
all.equal(distrib_cdf(d, 40, th, lower.tail = FALSE, log.p = TRUE),
          -(40 / 2)^1.5)
#> [1] TRUE

# There the survival probability itself has underflowed to zero.
distrib_cdf(d, 40, th, lower.tail = FALSE)
#> [1] 1.430608e-39
```
