# Gamma First Derivative in the Response, Mean and Dispersion

Computes the first derivative of the gamma log-density with respect to
the response, in closed form at the implied shape \\a = 1/\phi\\ and
rate \\b = 1/(\phi\mu)\\: \$\$\dfrac{\partial \ell}{\partial y} =
\dfrac{a - 1}{y} - b.\$\$ It changes sign at the mode \\y = (a-1)/b\\,
so it is positive below the mode and negative above it. At \\\phi = 1\\
the shape is 1, the first term drops out and the derivative is the
constant \\-1/\mu\\ of an exponential. At \\\phi \> 1\\ the shape falls
below 1, the density has no interior mode and the derivative is negative
throughout.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- y:

  A numeric vector of strictly positive observations. At `y = 0` the
  value is infinite unless the shape is exactly 1.

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

[`distrib_hess_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gamma1Distrib.md)
for the second derivative in the response,
[`distrib_gradient.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma1Distrib.md)
for the score in the parameters,
[`gamma1_shape_rate()`](https://statmodels7.github.io/distributions7/reference/gamma1_shape_rate.md)
for the conversion this uses, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)

# Written out at the implied shape and rate.
all.equal(distrib_grad_y(d, y, th), (2 - 1) / y - 1 / 1.5)
#> [1] TRUE

# Zero at the mode, positive below it and negative above.
mode <- (1 / 0.5 - 1) / (1 / (0.5 * 3))
c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#>    mode at_mode 
#>     1.5     0.0 
distrib_grad_y(d, c(0.5, 5), th)
#> [1]  1.3333333 -0.4666667

# At phi = 1 the family is exponential and the derivative is constant.
distrib_grad_y(d, y, list(mu = 3, phi = 1))
#> [1] -0.3333333 -0.3333333 -0.3333333

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  0.3333333 -0.3333333 -0.4666667
```
