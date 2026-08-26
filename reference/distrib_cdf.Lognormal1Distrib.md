# Lognormal Cumulative Distribution Function

Computes the lognormal distribution function \$\$F(q; \mu, \sigma^2) =
\Phi\left(\dfrac{\log q - \mu}{\sigma}\right), \qquad \sigma =
\sqrt{\sigma^2},\$\$ with \\\Phi\\ the standard normal distribution
function, by calling
[`stats::plnorm()`](https://rdrr.io/r/stats/Lognormal.html). The log
transformation is monotone, so this is the Gaussian distribution
function evaluated at \\\log q\\ and nothing else. Both tails are
available exactly, and `log.p = TRUE` returns a logarithm that stays
finite where the probability itself underflows.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below zero gives a
  lower-tail probability of 0.

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

[`distrib_quantile.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Lognormal1Distrib.md)
for the inverse,
[`distrib_pdf.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal1Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, which are the
Gaussian's at \\\log q\\ and so closed form, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
th <- list(mu = 0.5, sigma2 = 0.36)

# The method is stats::plnorm at meanlog = mu, sdlog = sqrt(sigma2).
all.equal(distrib_cdf(d, c(0.5, 1.6, 4), th),
          plnorm(c(0.5, 1.6, 4), 0.5, sqrt(0.36)))
#> [1] TRUE

# A monotone transformation, so it is the normal cdf at log q.
all.equal(distrib_cdf(d, c(0.5, 1.6, 4), th),
          pnorm(log(c(0.5, 1.6, 4)), 0.5, sqrt(0.36)))
#> [1] TRUE

# Half the mass lies below exp(mu), which is the median.
distrib_cdf(d, exp(0.5), th)
#> [1] 0.5

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 1e12, th, lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 1e12, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -1027.081
```
