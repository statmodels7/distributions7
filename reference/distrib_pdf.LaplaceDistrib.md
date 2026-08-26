# Laplace Probability Density Function

Computes the Laplace density \$\$f(y; \mu, \sigma) = \dfrac{1}{2\sigma}
\exp\left(-\dfrac{\|y - \mu\|}{\sigma}\right)\$\$ from the log-density
\\-\log(2\sigma) - \|y-\mu\|/\sigma\\, which is formed first and
exponentiated only when `log = FALSE`. The density is continuous
everywhere and has a corner at \\y = \mu\\, where the two exponential
arms meet. That corner makes the family non-regular in its location.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive; the arithmetic is
  performed as written, so a non-positive value gives `NaN` or `Inf`
  without a warning of its own.

- log:

  Logical of length 1. When `TRUE` the log-density is returned, which is
  exact and finite at every real `y`. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LaplaceDistrib.md)
for the distribution function,
[`distrib_grad_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LaplaceDistrib.md)
for the kink at \\y = \mu\\,
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
for the rate parametrization, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

# The density, written out.
all.equal(distrib_pdf(d, y, th), exp(-abs(y - 0.4) / 1.5) / (2 * 1.5))
#> [1] TRUE

# The log-density is exactly linear in |y - mu|, so it is a pair of
# straight lines meeting at mu.
round(distrib_pdf(d, 0.4 + c(-3, -1.5, 0, 1.5, 3), th, log = TRUE), 6)
#> [1] -3.098612 -2.098612 -1.098612 -2.098612 -3.098612

# Symmetric about mu, and at its maximum there.
distrib_pdf(d, 0.4 + c(-2, 2), th)
#> [1] 0.08786571 0.08786571
c(at_mu = distrib_pdf(d, 0.4, th), one_over_2sigma = 1 / (2 * 1.5))
#>           at_mu one_over_2sigma 
#>       0.3333333       0.3333333 
```
