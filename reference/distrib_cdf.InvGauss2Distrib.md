# Inverse Gaussian Cumulative Distribution Function in Mean and Shape

Computes the inverse Gaussian distribution function, which is elementary
in the standard normal distribution function \\\Phi\\: \$\$F(q; \mu,
\lambda) = \Phi\left\\\sqrt{\dfrac{\lambda}{q}} \left(\dfrac{q}{\mu} -
1\right)\right\\ + e^{2\lambda/\mu}\\
\Phi\left\\-\sqrt{\dfrac{\lambda}{q}} \left(\dfrac{q}{\mu} +
1\right)\right\\,\$\$ by calling
[`statmod::pinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = 1/lambda`. Both tails are available
exactly, and `log.p = TRUE` returns a logarithm that stays finite where
the probability itself underflows. The exponential factor overflows at
ordinary settings, so that function evaluates the expression on the log
scale.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below zero gives a
  lower-tail probability of 0.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
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
`max(length(q), length(mu), length(lambda))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGauss2Distrib.md)
for the inverse,
[`distrib_pdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGauss2Distrib.md)
for the density,
[`distrib_grad_cdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.InvGauss2Distrib.md)
for the derivatives of this function in the parameters, which are closed
form here, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
th <- list(mu = 2, lambda = 3)

# The method is statmod::pinvgauss at dispersion 1/lambda.
all.equal(distrib_cdf(d, c(1, 2, 3), th),
          statmod::pinvgauss(c(1, 2, 3), mean = 2, dispersion = 1 / 3))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, 2, th) + distrib_cdf(d, 2, th, lower.tail = FALSE)
#> [1] 1

# The law is right skewed, so most of the mass sits below the mean.
distrib_cdf(d, 2, th)
#> [1] 0.6436706

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 1e4, th, lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 1e4, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -3761.705
```
