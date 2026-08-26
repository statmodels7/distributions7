# Laplace Probability Density Function, Rate Parametrization

Computes the Laplace density in the rate parametrization, \$\$f(y; \mu,
\lambda) = \dfrac{\lambda}{2} \exp\left(-\lambda\|y - \mu\|\right),\$\$
from the log-density \\\log(\lambda/2) - \lambda\|y-\mu\|\\, which is
formed first and exponentiated only when `log = FALSE`. It equals the
density of
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
at \\\sigma = 1/\lambda\\.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `lambda` must be strictly positive; the arithmetic is
  performed as written, so a non-positive value gives `NaN` without a
  warning of its own.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(lambda))`, one value per observation.

## See also

[`distrib_pdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.LaplaceDistrib.md)
for the scale parametrization,
[`distrib_cdf.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Laplace2Distrib.md)
for the distribution function, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)

# The density, written out.
all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 2)),
          exp(-2 * abs(y - 0.4)) * 2 / 2)
#> [1] TRUE

# The same law as the scale parametrization at lambda = 1/sigma.
all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 1 / 1.5)),
          distrib_pdf(laplace_distrib(), y, list(mu = 0.4, sigma = 1.5)))
#> [1] TRUE

# A larger rate concentrates the mass: the peak is lambda/2.
vapply(c(0.5, 1, 2), function(l) distrib_pdf(d, 0.4,
                                             list(mu = 0.4, lambda = l)),
       numeric(1))
#> [1] 0.25 0.50 1.00
```
