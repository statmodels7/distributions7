# Gaussian Cumulative Distribution Function

Computes the Gaussian distribution function \$\$F(q; \mu, \sigma) =
\Phi\left(\dfrac{q-\mu}{\sigma}\right)\$\$ with \\\Phi\\ the standard
normal distribution function, by calling
[`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html). Both tails are
available exactly: `lower.tail = FALSE` evaluates \\1 - F\\ without
forming the difference, and `log.p = TRUE` returns a logarithm that
stays finite where the probability itself underflows to zero.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

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

[`distrib_quantile.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gaussian1Distrib.md)
for the inverse,
[`distrib_pdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian1Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- gaussian1_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The method is stats::pnorm at this parametrization.
all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
          pnorm(c(-1.2, 0.3, 2.5), mean = 0.4, sd = 1.5))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
#> [1] 1

# Forty standard deviations out the upper tail underflows; its log does not.
distrib_cdf(d, 40, list(mu = 0, sigma = 1), lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 40, list(mu = 0, sigma = 1), lower.tail = FALSE, log.p = TRUE)
#> [1] -804.6084
```
