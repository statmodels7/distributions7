# Inverse Gaussian Cumulative Distribution Function in Mean and Dispersion

Computes the inverse Gaussian distribution function, which is elementary
in the standard normal distribution function \\\Phi\\: \$\$F(q; \mu,
\phi) = \Phi\left\\\sqrt{\dfrac{1}{\phi q}} \left(\dfrac{q}{\mu} -
1\right)\right\\ + e^{2/(\phi\mu)}\\ \Phi\left\\-\sqrt{\dfrac{1}{\phi
q}} \left(\dfrac{q}{\mu} + 1\right)\right\\,\$\$ by calling
[`statmod::pinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = phi`. Both tails are available exactly,
and `log.p = TRUE` returns a logarithm that stays finite where the
probability itself underflows.

The exponential factor is the reason this expression is evaluated on the
log scale rather than as written: \\e^{2/(\phi\mu)}\\ overflows at
ordinary settings, at \\\mu = 0.01\\ and \\\phi = 0.1\\ the exponent
already being 2000, and it multiplies a normal tail that underflows by
the same amount.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- q:

  A numeric vector of quantiles. A value at or below zero gives a
  lower-tail probability of 0.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `q`. A component of length 1 is recycled.
  Both must be strictly positive.

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
`max(length(q), length(mu), length(phi))`. With `log.p = TRUE` the
values are logarithms and are non-positive.

## See also

[`distrib_quantile.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.InvGauss1Distrib.md)
for the inverse,
[`distrib_pdf.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGauss1Distrib.md)
for the density,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the derivatives of this function in the parameters, which are closed
form here because the expression above is elementary, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
th <- list(mu = 1, phi = 2)
q <- c(0.5, 1, 2)

# The method is statmod::pinvgauss at this parametrization.
all.equal(distrib_cdf(d, q, th),
          statmod::pinvgauss(q, mean = 1, dispersion = 2))
#> [1] TRUE

# The closed form above, evaluated directly at these safe values.
a <- sqrt(1 / (2 * q)) * (q / 1 - 1)
b <- -sqrt(1 / (2 * q)) * (q / 1 + 1)
pnorm(a) + exp(2 / (2 * 1)) * pnorm(b)
#> [1] 0.4901383 0.7137918 0.8730633

# The two tails sum to one.
distrib_cdf(d, 2, th) + distrib_cdf(d, 2, th, lower.tail = FALSE)
#> [1] 1

# The law is heavily right skewed, so most of the mass sits below the mean.
distrib_cdf(d, 1, th)
#> [1] 0.7137918

# Far in the upper tail the probability underflows and its log does not.
distrib_cdf(d, 1e4, th, lower.tail = FALSE)
#> [1] 0
distrib_cdf(d, 1e4, th, lower.tail = FALSE, log.p = TRUE)
#> [1] -2513.195
```
