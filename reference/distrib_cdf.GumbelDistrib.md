# Gumbel Cumulative Distribution Function

Computes the Gumbel distribution function, with \\z = (q -
\mu)/\sigma\\, \$\$F(q; \mu, \sigma) = \exp\left\\-e^{-z}\right\\.\$\$
This is the defining property of the family: a maximum of \\n\\
independent light-tailed variables has, after centering and scaling, a
distribution function that is the \\n\\th power of one distribution
function, and \\\exp\\-e^{-z}\\\\ is the fixed point of that operation.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

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

## Details

Both tails are computed on the scale that keeps them accurate. The lower
tail is exact on the log scale, \\\log F = -e^{-z}\\, and the upper tail
goes through [`base::expm1()`](https://rdrr.io/r/base/Log.html), so
neither loses precision where it is small.

## See also

[`distrib_quantile.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GumbelDistrib.md)
for the inverse,
[`distrib_pdf.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GumbelDistrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
th <- list(mu = 0, sigma = 1)
q <- c(-1, 0, 1)

# The closed form, written out.
all.equal(distrib_cdf(d, q, th), exp(-exp(-q)))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 1, th) + distrib_cdf(d, 1, th, lower.tail = FALSE)
#> [1] 1

# Max-stability: the cdf raised to the nth power is the same law shifted
# by sigma log n.
n <- 10
all.equal(distrib_cdf(d, q, th)^n,
          distrib_cdf(d, q, list(mu = log(n), sigma = 1)))
#> [1] TRUE

# Deep in the lower tail the probability underflows and its log does not.
distrib_cdf(d, -700, th)
#> [1] 0
distrib_cdf(d, -700, th, log.p = TRUE)
#> [1] -1.014232e+304
```
