# Inverse Gaussian Probability Density Function in Mean and Shape

Computes the inverse Gaussian density \$\$f(y; \mu, \lambda) =
\sqrt{\dfrac{\lambda}{2\pi y^3}}
\exp\left\\-\dfrac{\lambda(y-\mu)^2}{2\mu^2 y}\right\\, \qquad y \>
0,\$\$ by calling
[`statmod::dinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = 1/lambda`, that package writing the
family in the dispersion. With `log = TRUE` the logarithm is formed
inside `dinvgauss()` and stays finite where the density itself
underflows.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, \infty)\\; a
  value at or below zero gives 0.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(lambda))`, one value per observation.

## See also

[`distrib_cdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGauss2Distrib.md)
for the distribution function,
[`distrib_gradient.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss2Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.InvGauss1Distrib.md)
for the same density in the dispersion, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
y <- c(1, 2, 3)
th <- list(mu = 2, lambda = 3)

# The method is statmod::dinvgauss at dispersion 1/lambda.
all.equal(distrib_pdf(d, y, th),
          statmod::dinvgauss(y, mean = 2, dispersion = 1 / 3))
#> [1] TRUE

# The same law as invgauss1 at phi = 1/lambda.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(invgauss1_distrib(), y, list(mu = 2, phi = 1 / 3)))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(1, 2, 4), lambda = 3))
#> [1] 0.6909883 0.2443013 0.1288894

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, 1e4, th)
#> [1] 0
distrib_pdf(d, 1e4, th, log = TRUE)
#> [1] -3762.685
```
