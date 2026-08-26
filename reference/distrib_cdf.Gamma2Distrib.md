# Gamma Cumulative Distribution Function in Mean and Variance

Computes the gamma distribution function, the regularized incomplete
gamma function \$\$F(q; \mu, \sigma^2) = \dfrac{\gamma(\alpha, \lambda
q)}{\Gamma(\alpha)}, \qquad \alpha = \dfrac{\mu^2}{\sigma^2}, \quad
\lambda = \dfrac{\mu}{\sigma^2},\$\$ with \\\gamma\\ the lower
incomplete gamma function, by calling
[`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html) at that
shape and rate. Both tails are available exactly: `lower.tail = FALSE`
evaluates \\1 - F\\ without forming the difference, and `log.p = TRUE`
returns a logarithm that stays finite where the probability itself
underflows to zero.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below zero gives a
  lower-tail probability of 0.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `q`. A component of length 1 is
  recycled. Both must be strictly positive.

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

[`distrib_quantile.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.Gamma2Distrib.md)
for the inverse,
[`distrib_pdf.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma2Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, which the gamma
takes by finite difference because the derivative of an incomplete gamma
in its shape is hypergeometric, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
th <- list(mu = 3, sigma2 = 2)

# The method is stats::pgamma at the implied shape and rate.
all.equal(distrib_cdf(d, c(1, 3, 5), th),
          pgamma(c(1, 3, 5), shape = 9 / 2, rate = 3 / 2))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 5, th) + distrib_cdf(d, 5, th, lower.tail = FALSE)
#> [1] 1

# The gamma is right skewed, so less than half the mass lies below the mean.
distrib_cdf(d, 3, th)
#> [1] 0.5627258

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 3000, th, lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 3000, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -4473.012
```
