# Cauchy Probability Density Function

Computes the Cauchy density \$\$f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma
\left\[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right\]}\$\$ by calling
[`stats::dcauchy()`](https://rdrr.io/r/stats/Cauchy.html) at
`location = mu` and `scale = sigma`. The density decays like \\y^{-2}\\,
so it stays representable far out where a Gaussian of the same scale has
underflowed to zero.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive; a zero or negative value
  gives `NaN` with a warning from
  [`stats::dcauchy()`](https://rdrr.io/r/stats/Cauchy.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.CauchyDistrib.md)
for the distribution function,
[`distrib_gradient.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian1Distrib.md)
for the light-tailed comparison, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
y <- c(-1.2, 0.3, 2.5)

# The method is stats::dcauchy at this parametrization.
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)),
          dcauchy(y, location = 0.4, scale = 1.5))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0, 1, 2), sigma = c(1, 1.5, 2)))
#> [1] 0.1304549 0.1742572 0.1497929

# Forty scale units out the Cauchy density is still an ordinary number,
# where the Gaussian of the same scale has underflowed.
distrib_pdf(d, 40, list(mu = 0, sigma = 1))
#> [1] 0.0001988194
distrib_pdf(gaussian1_distrib(), 40, list(mu = 0, sigma = 1))
#> [1] 0
```
