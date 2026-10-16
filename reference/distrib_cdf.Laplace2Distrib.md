# Laplace Cumulative Distribution Function, Rate Parametrization

Computes the Laplace distribution function in the rate parametrization,
one exponential on each side of the location: \$\$F(q; \mu, \lambda) =
\begin{cases} \tfrac{1}{2}\exp(\lambda(q-\mu)), & q \< \mu,\\\[4pt\] 1 -
\tfrac{1}{2}\exp(-\lambda(q-\mu)), & q \ge \mu. \end{cases}\$\$ Both
arms are continuous at \\q = \mu\\, where the value is \\1/2\\. Each
tail is computed on the log scale at \\r = \lambda(q - \mu)\\ by
[`laplace_cdf_tails()`](https://statmodels7.github.io/distributions7/reference/laplace_cdf_tails.md),
as for
[`distrib_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LaplaceDistrib.md):
the far tail is \\\log(1/2) - \|r\|\\ and the near one \\\log(1 -
e^{-\|r\|}/2)\\, so `lower.tail = FALSE` and `log.p = TRUE` subtract
from nothing, and a tail probability far below machine epsilon, or its
logarithm, is returned accurately.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. `lambda` must be strictly positive.

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
`max(length(q), length(mu), length(lambda))`.

## See also

[`distrib_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LaplaceDistrib.md)
for the scale parametrization,
[`distrib_quantile.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Laplace2Distrib.md)
for the inverse, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
th <- list(mu = 0.4, lambda = 2)

# The distribution function equals one half at the location.
distrib_cdf(d, 0.4, th)
#> [1] 0.5

# Each arm is an exponential, written out.
q <- c(-2, 3)
all.equal(distrib_cdf(d, q, th),
          ifelse(q < 0.4, 0.5 * exp(2 * (q - 0.4)),
                 1 - 0.5 * exp(-2 * (q - 0.4))))
#> [1] TRUE

# The probabilities equal those of the scale parametrization at
# sigma = 1/lambda.
all.equal(distrib_cdf(d, q, list(mu = 0.4, lambda = 1 / 1.5)),
          distrib_cdf(laplace_distrib(), q, list(mu = 0.4, sigma = 1.5)))
#> [1] TRUE
```
