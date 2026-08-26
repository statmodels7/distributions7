# Chi-Squared Cumulative Distribution Function

Computes the chi-squared distribution function, the regularized
incomplete gamma function at shape \\\mu/2\\ and argument \\q/2\\, by
calling [`stats::pchisq()`](https://rdrr.io/r/stats/Chisquare.html) at
`df = mu`. Both tails are available exactly: `lower.tail = FALSE`
evaluates \\1 - F\\ without forming the difference, and `log.p = TRUE`
returns a logarithm that stays finite where the probability itself
underflows to zero.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below zero gives a
  lower-tail probability of 0.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `q`, recycled if of length 1. It must be strictly
  positive.

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
`max(length(q), length(mu))`. With `log.p = TRUE` the values are
logarithms and are non-positive.

## See also

[`distrib_quantile.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ChisqDistrib.md)
for the inverse,
[`distrib_pdf.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ChisqDistrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameter, which the
chi-squared takes by finite difference because the derivative of an
incomplete gamma in its shape is hypergeometric, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
th <- list(mu = 4)

# The method is stats::pchisq at df = mu.
all.equal(distrib_cdf(d, c(1, 4, 9), th), pchisq(c(1, 4, 9), df = 4))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 9, th) + distrib_cdf(d, 9, th, lower.tail = FALSE)
#> [1] 1

# Right skewed, so less than half the mass lies below the mean.
distrib_cdf(d, 4, th)
#> [1] 0.5939942

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 2000, th, lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 2000, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -993.0912
```
