# Laplace Cumulative Distribution Function

Computes the Laplace distribution function, which is one exponential on
each side of the location: \$\$F(q; \mu, \sigma) = \begin{cases}
\tfrac{1}{2}\exp\\\left(\dfrac{q-\mu}{\sigma}\right), & q \<
\mu,\\\[4pt\] 1 - \tfrac{1}{2}\exp\\\left(-\dfrac{q-\mu}{\sigma}\right),
& q \ge \mu. \end{cases}\$\$ The branch is selected per observation, and
both arms are continuous at \\q = \mu\\, where the value is \\1/2\\.
With `lower.tail = FALSE` the complement is taken after the branch, and
with `log.p = TRUE` the logarithm of the result is taken; neither is
computed in a way that avoids cancellation, so a probability that has
already underflowed is not recovered.

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
  \le q)\\; when `FALSE` they are \\P(Y \> q)\\, computed as \\1 - F\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned, taken after the probability itself is formed. Defaults to
  `FALSE`.

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

# One half at the location, the two arms meeting there.
distrib_cdf(d, 0.4, th)
#> [1] 0.5

# Each arm is an exponential, written out.
q <- c(-2, 3)
all.equal(distrib_cdf(d, q, th),
          ifelse(q < 0.4, 0.5 * exp((q - 0.4) / 1.5),
                 1 - 0.5 * exp(-(q - 0.4) / 1.5)))
#> [1] TRUE

# Symmetric: F(mu - a) + F(mu + a) = 1.
distrib_cdf(d, 0.4 - 2, th) + distrib_cdf(d, 0.4 + 2, th)
#> [1] 1
```
