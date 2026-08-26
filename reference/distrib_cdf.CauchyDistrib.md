# Cauchy Cumulative Distribution Function

Computes the Cauchy distribution function \$\$F(q; \mu, \sigma) =
\dfrac{1}{2} +
\dfrac{1}{\pi}\arctan\left(\dfrac{q-\mu}{\sigma}\right)\$\$ by calling
[`stats::pcauchy()`](https://rdrr.io/r/stats/Cauchy.html). Far out in
either tail the probability behaves like \\1/(\pi \|z\|)\\ with \\z =
(q-\mu)/\sigma\\, so it decays at a polynomial rate and
`lower.tail = FALSE` returns an ordinary number where a light-tailed
family returns zero.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

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

[`distrib_quantile.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.CauchyDistrib.md)
for the inverse,
[`distrib_pdf.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.CauchyDistrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The method is stats::pcauchy at this parametrization.
all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
          pcauchy(c(-1.2, 0.3, 2.5), location = 0.4, scale = 1.5))
#> [1] TRUE

# A hundred scale units out the tail is 1/(100 pi), not zero.
distrib_cdf(d, 0.4 + 100 * 1.5, th, lower.tail = FALSE)
#> [1] 0.003182993
1 / (100 * pi)
#> [1] 0.003183099

# The two tails sum to one.
distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
#> [1] 1
```
