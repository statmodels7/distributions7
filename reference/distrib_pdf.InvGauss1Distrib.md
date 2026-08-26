# Inverse Gaussian Probability Density Function in Mean and Dispersion

Computes the inverse Gaussian density \$\$f(y; \mu, \phi) =
\sqrt{\dfrac{1}{2\pi\phi y^3}} \exp\left\\-\dfrac{(y-\mu)^2}{2\phi\mu^2
y}\right\\, \qquad y \> 0,\$\$ by calling
[`statmod::dinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = phi`. With `log = TRUE` the logarithm
is formed inside that function and stays finite where the density itself
underflows, which happens quickly: the density falls off like
\\\exp\\-y/(2\phi\mu^2)\\\\ in the right tail and like \\\exp\\-1/(2\phi
y)\\\\ at the origin.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, \infty)\\; a
  value at or below zero gives 0.

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

[`distrib_cdf.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGauss1Distrib.md)
for the distribution function,
[`distrib_gradient.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss1Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGauss2Distrib.md)
for the same density in the variance, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)

# The method is statmod::dinvgauss at this parametrization.
all.equal(distrib_pdf(d, y, th),
          statmod::dinvgauss(y, mean = 1, dispersion = 2))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0.5, 1, 2), phi = 2))
#> [1] 0.79788456 0.28209479 0.09973557

# The density vanishes at the origin faster than any power.
distrib_pdf(d, c(1e-2, 1e-3, 1e-4), th)
#> [1]  6.443095e-09 3.924761e-105  0.000000e+00

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, 1e4, th)
#> [1] 0
distrib_pdf(d, 1e4, th, log = TRUE)
#> [1] -2514.581
```
