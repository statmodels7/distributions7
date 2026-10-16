# Laplace Cumulative Distribution Function

Computes the Laplace distribution function, which is one exponential on
each side of the location: \$\$F(q; \mu, \sigma) = \begin{cases}
\tfrac{1}{2}\exp\\\left(\dfrac{q-\mu}{\sigma}\right), & q \<
\mu,\\\[4pt\] 1 - \tfrac{1}{2}\exp\\\left(-\dfrac{q-\mu}{\sigma}\right),
& q \ge \mu. \end{cases}\$\$ Both arms are continuous at \\q = \mu\\,
where the value is \\1/2\\. The method computes the logarithm of each
tail at the standardized residual \\r = (q-\mu)/\sigma\\: the tail on
the far side of the location is \\\log(1/2) - \|r\|\\, and the tail on
the near side is \\\log(1 - e^{-\|r\|}/2)\\, computed with
[`log1p()`](https://rdrr.io/r/base/Log.html). The requested tail is
returned on the log scale or exponentiated, so the method does not
subtract from one when `lower.tail = FALSE` or `log.p = TRUE`, and a
tail probability far below machine epsilon, or its logarithm, is
returned accurately. The helper is
[`laplace_cdf_tails()`](https://statmodels7.github.io/distributions7/reference/laplace_cdf_tails.md).

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

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

[`distrib_quantile.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LaplaceDistrib.md)
for the inverse,
[`distrib_pdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LaplaceDistrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The value is one half at the location, where the two arms meet.
distrib_cdf(d, 0.4, th)
#> [1] 0.5

# Each arm is an exponential, written out.
q <- c(-2, 3)
all.equal(distrib_cdf(d, q, th),
          ifelse(q < 0.4, 0.5 * exp((q - 0.4) / 1.5),
                 1 - 0.5 * exp(-(q - 0.4) / 1.5)))
#> [1] TRUE

# The distribution is symmetric: F(mu - a) + F(mu + a) = 1.
distrib_cdf(d, 0.4 - 2, th) + distrib_cdf(d, 0.4 + 2, th)
#> [1] 1

# Fifty scale units out, the upper tail is exp(-50)/2, and its logarithm
# is exact as well.
distrib_cdf(d, 0.4 + 50 * 1.5, th, lower.tail = FALSE)
#> [1] 9.643749e-23
distrib_cdf(d, 0.4 + 50 * 1.5, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -50.69315
```
