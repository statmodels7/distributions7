# Inverse Gaussian First Derivative in the Response

Computes the first derivative of the inverse Gaussian log-density with
respect to the response, in closed form: \$\$\dfrac{\partial
\ell}{\partial y} = -\dfrac{3}{2y} - \dfrac{y^2 - \mu^2}{2\phi\mu^2
y^2}.\$\$ The first term comes from the \\y^{-3/2}\\ factor of the
density and the second from the exponent. It changes sign at the mode,
which for this family is \\\mu\\(1 + 9\phi^2\mu^2/4)^{1/2} -
3\phi\mu/2\\\\, always strictly below the mean.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of strictly positive observations. The value diverges
  as `y` approaches zero.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(phi))`,
one value per observation.

## See also

[`distrib_hess_y.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.InvGauss1Distrib.md)
for the second derivative in the response,
[`distrib_gradient.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss1Distrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)

# Written out.
all.equal(distrib_grad_y(d, y, th),
          -1.5 / y - (y^2 - 1^2) / (2 * 2 * 1^2 * y^2))
#> [1] TRUE

# Zero at the mode, which lies below the mean.
mode <- 1 * (sqrt(1 + 9 * 2^2 * 1^2 / 4) - 3 * 2 * 1 / 2)
c(mode = mode, mean = 1, at_mode = distrib_grad_y(d, mode, th))
#>          mode          mean       at_mode 
#>  1.622777e-01  1.000000e+00 -1.243450e-14 

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1] -2.2500 -1.5000 -0.9375
```
