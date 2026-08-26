# Logistic Cumulative Distribution Function

Computes the logistic distribution function \$\$F(q; \mu, \sigma) =
\dfrac{1}{1 + \exp\left(-\dfrac{q-\mu}{\sigma}\right)}\$\$ by calling
[`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html). This is the
logistic sigmoid, the inverse of the logit link, so the same curve
appears here as a distribution function and in `linkfunctions7` as an
inverse link.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

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
  returned, which stays finite far into either tail. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(q), length(mu), length(sigma))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.LogisticDistrib.md)
for the inverse,
[`distrib_pdf.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LogisticDistrib.md)
for the density,
[`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
whose inverse is this curve, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The method is stats::plogis at this parametrization.
all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
          plogis(c(-1.2, 0.3, 2.5), location = 0.4, scale = 1.5))
#> [1] TRUE

# Symmetric: F(mu - a) + F(mu + a) = 1.
distrib_cdf(d, 0.4 - 2, th) + distrib_cdf(d, 0.4 + 2, th)
#> [1] 1

# It is the inverse logit link, so the two agree exactly.
all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), list(mu = 0, sigma = 1)),
          linkfunctions7::linkinv(linkfunctions7::logit_link(),
                                  c(-1.2, 0.3, 2.5)))
#> [1] TRUE
```
