# Gaussian Cumulative Distribution Function in Mean and Variance

Computes the Gaussian distribution function \$\$F(q; \mu, \sigma^2) =
\Phi\left(\dfrac{q-\mu}{\sigma}\right), \qquad \sigma =
\sqrt{\sigma^2},\$\$ with \\\Phi\\ the standard normal distribution
function, by calling
[`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html). Both tails are
available exactly: `lower.tail = FALSE` evaluates \\1 - F\\ without
forming the difference, and `log.p = TRUE` returns a logarithm that
stays finite where the probability itself underflows to zero.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

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
`max(length(q), length(mu), length(sigma2))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gaussian2Distrib.md)
for the inverse,
[`distrib_pdf.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian2Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- gaussian2_distrib()
th <- list(mu = 1, sigma2 = 4)

# The method is stats::pnorm at sd = sqrt(sigma2).
all.equal(distrib_cdf(d, c(-1.2, 0.3, 2.5), th),
          pnorm(c(-1.2, 0.3, 2.5), mean = 1, sd = 2))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 3, th) + distrib_cdf(d, 3, th, lower.tail = FALSE)
#> [1] 1

# Forty standard deviations out the upper tail underflows; its log does not.
distrib_cdf(d, 40, list(mu = 0, sigma2 = 1), lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 40, list(mu = 0, sigma2 = 1), lower.tail = FALSE,
            log.p = TRUE)
#> [1] -804.6084
```
