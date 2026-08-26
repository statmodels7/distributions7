# Exponential Cumulative Distribution Function

Computes the exponential distribution function \$\$F(q; \mu) = 1 -
\exp\left(-\dfrac{q}{\mu}\right), \qquad q \ge 0,\$\$ by calling
[`stats::pexp()`](https://rdrr.io/r/stats/Exponential.html) at
`rate = 1/mu`. With `lower.tail = FALSE` the survival function
\\\exp(-q/\mu)\\ is returned directly, without forming the difference,
so it stays exact far into the tail; combined with `log.p = TRUE` it is
simply \\-q/\mu\\.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- q:

  A numeric vector of quantiles. A negative value gives 0.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `q`. `mu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\; when `FALSE` they are the survival function \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu))`.

## See also

[`distrib_quantile.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.ExponentialDistrib.md)
for the inverse,
[`distrib_pdf.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ExponentialDistrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
th <- list(mu = 2)

# The method is stats::pexp at rate = 1/mu.
all.equal(distrib_cdf(d, c(0.3, 1.1, 4.0), th),
          pexp(c(0.3, 1.1, 4.0), rate = 1 / 2))
#> [1] TRUE

# The survival function on the log scale is exactly -q/mu.
distrib_cdf(d, c(10, 100, 1000), th, lower.tail = FALSE, log.p = TRUE)
#> [1]   -5  -50 -500
-c(10, 100, 1000) / 2
#> [1]   -5  -50 -500

# Memoryless: the chance of surviving another unit does not depend on how
# long the wait has already been.
c(distrib_cdf(d, 4, th, lower.tail = FALSE) /
    distrib_cdf(d, 3, th, lower.tail = FALSE),
  distrib_cdf(d, 1, th, lower.tail = FALSE))
#> [1] 0.6065307 0.6065307
```
