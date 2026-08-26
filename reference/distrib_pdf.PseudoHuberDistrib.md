# Pseudo-Huber Probability Density Function

Computes the pseudo-Huber density \$\$f(y; \mu, \sigma, \nu) =
\dfrac{1}{2 \sigma \sqrt{\nu}\\ K_1(\sqrt{\nu})} \exp\left(-\sqrt{\nu +
\left(\dfrac{y-\mu}{\sigma}\right)^2}\right),\$\$ with \\K_1\\ the
modified Bessel function of the second kind. The exponent is the
negative pseudo-Huber loss, quadratic in the residual near the location
and linear far from it. That exponent is the pseudo-Huber loss, and this
density is its exponential.

The normalizing constant is formed on the log scale through the
**exponentially scaled** Bessel function, \\\log K_1(x) = \log\\e^x
K_1(x)\\ - x\\. The Bessel terms are degree-homogeneous, so the scaled
form is exact and stays finite where \\K_1(\sqrt{\nu})\\ itself
underflows, verified to \\\nu = 2000\\.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma), length(nu))`, one value per
observation.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
shape and \\K_1\\ the modified Bessel function of the second kind of
order one, `besselK(x, 1)` in R.

## See also

[`distrib_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuberDistrib.md)
for the distribution function,
[`distrib_gradient.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md)
for the derivatives of the log-density,
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the two limits, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The formula written out, with the Bessel constant.
D <- sqrt(2 + ((y - 0.4) / 1.2)^2)
all.equal(distrib_pdf(d, y, th),
          exp(-D) / (2 * 1.2 * sqrt(2) * besselK(sqrt(2), 1)))
#> [1] TRUE

# It integrates to one.
integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#> [1] 1

# A small shape is a Laplace of scale sigma; a large one a Gaussian of
# standard deviation sigma * nu^(1/4).
yy <- c(0.5, 1, 2, 4)
max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e-8)) -
        0.5 * exp(-abs(yy))))
#> [1] 1.186719e-08
max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e8)) -
        dnorm(yy, 0, 100)))
#> [1] 1.495912e-07
```
