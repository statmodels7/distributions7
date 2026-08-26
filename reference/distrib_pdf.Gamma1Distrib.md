# Gamma Probability Density Function in Mean and Dispersion

Computes the gamma density \$\$f(y; \mu, \phi) = \dfrac{y^{1/\phi -
1}\\e^{-y/(\phi\mu)}} {(\phi\mu)^{1/\phi}\\\Gamma(1/\phi)}, \qquad y \>
0,\$\$ by calling
[`stats::dgamma()`](https://rdrr.io/r/stats/GammaDist.html) at shape \\a
= 1/\phi\\ and rate \\b = 1/(\phi\mu)\\. With `log = TRUE` the logarithm
is formed inside [`dgamma()`](https://rdrr.io/r/stats/GammaDist.html)
and stays finite where the density itself underflows.

The density is unbounded at the origin when \\\phi \> 1\\, where the
shape falls below one; it is flat there at \\\phi = 1\\, the exponential
case, and vanishes at the origin for \\\phi \< 1\\.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, \infty)\\; a
  negative value gives 0 and `y = 0` gives 0, `Inf` or the rate
  according to whether \\\phi\\ is below, above or equal to 1.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(phi))`, one value per observation.

## See also

[`distrib_cdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma1Distrib.md)
for the distribution function,
[`distrib_gradient.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma1Distrib.md)
for the derivatives of the log-density,
[`gamma1_shape_rate()`](https://statmodels7.github.io/distributions7/reference/gamma1_shape_rate.md)
for the conversion this uses, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)

# The method is stats::dgamma at shape 1/phi and rate 1/(phi mu).
all.equal(distrib_pdf(d, y, th),
          dgamma(y, shape = 1 / 0.5, rate = 1 / (0.5 * 3)))
#> [1] TRUE

# At phi = 1 the shape is 1 and the gamma is the exponential.
all.equal(distrib_pdf(d, y, list(mu = 3, phi = 1)),
          distrib_pdf(exponential_distrib(), y, list(mu = 3)))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(1, 3, 9), phi = 0.5))
#> [1] 0.54134113 0.18044704 0.08128222

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, 1e4, th)
#> [1] 0
distrib_pdf(d, 1e4, th, log = TRUE)
#> [1] -6658.267
```
