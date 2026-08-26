# Weibull Probability Density Function

Computes the Weibull density \$\$f(y; \mu, \sigma) = \dfrac{\sigma}{\mu}
\left(\dfrac{y}{\mu}\right)^{\sigma - 1}
\exp\left\\-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\\, \qquad y \>
0,\$\$ by calling
[`stats::dweibull()`](https://rdrr.io/r/stats/Weibull.html) at
`shape = sigma` and `scale = mu`, so the accuracy and the underflow
behavior are R's own. Outside the support the density is 0. With
`log = TRUE` the logarithm is formed inside
[`dweibull()`](https://rdrr.io/r/stats/Weibull.html) and stays finite
far into the upper tail, where the density itself underflows.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of observations. A value at or below zero is outside
  the support and gives a density of 0.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive; a zero or negative value
  gives `NaN` with a warning from
  [`stats::dweibull()`](https://rdrr.io/r/stats/Weibull.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Weibull1Distrib.md)
for the distribution function,
[`distrib_gradient.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull1Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic and
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
for the family.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)

# The method is stats::dweibull with the shape and the scale swapped into
# this parametrization's order.
all.equal(distrib_pdf(d, y, list(mu = 2, sigma = 1.5)),
          dweibull(y, shape = 1.5, scale = 2))
#> [1] TRUE

# Shape 1 is the exponential with mean mu, and shape 2 the Rayleigh.
all.equal(distrib_pdf(d, y, list(mu = 2, sigma = 1)), dexp(y, rate = 1 / 2))
#> [1] TRUE

# Below the support the density is zero; in the far upper tail it
# underflows and its logarithm does not.
distrib_pdf(d, c(-1, 0), list(mu = 2, sigma = 1.5))
#> [1] 0 0
distrib_pdf(d, 60, list(mu = 2, sigma = 1.5))
#> [1] 1.785487e-71
distrib_pdf(d, 60, list(mu = 2, sigma = 1.5), log = TRUE)
#> [1] -162.9039
```
