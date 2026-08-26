# Gamma First Derivative in the Response, Mean and Variance

Computes the first derivative of the gamma log-density with respect to
the response, in closed form at the implied shape \\\alpha =
\mu^2/\sigma^2\\ and rate \\\lambda = \mu/\sigma^2\\:
\$\$\dfrac{\partial \ell}{\partial y} = \dfrac{\alpha - 1}{y} -
\lambda.\$\$ It changes sign at the mode \\y = (\alpha-1)/\lambda\\, so
it is positive below the mode and negative above it. When \\\sigma^2 =
\mu^2\\ the shape is 1, the first term drops out and the derivative is
the constant \\-1/\mu\\ of an exponential; when \\\sigma^2 \> \mu^2\\
the density has no interior mode and the derivative is negative
throughout.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of strictly positive observations. At `y = 0` the
  value is infinite unless the shape is exactly 1.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma2))`,
one value per observation.

## See also

[`distrib_hess_y.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gamma2Distrib.md)
for the second derivative in the response,
[`distrib_gradient.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma2Distrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)

# Written out at the implied shape and rate.
all.equal(distrib_grad_y(d, y, th), (9 / 2 - 1) / y - 3 / 2)
#> [1] TRUE

# Zero at the mode, positive below it and negative above.
mode <- (9 / 2 - 1) / (3 / 2)
c(mode = mode, at_mode = distrib_grad_y(d, mode, th))
#>     mode  at_mode 
#> 2.333333 0.000000 

# At sigma2 = mu^2 the family is exponential and the derivative is constant.
distrib_grad_y(d, y, list(mu = 3, sigma2 = 9))
#> [1] -0.3333333 -0.3333333 -0.3333333

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  2.0000000 -0.3333333 -0.8000000
```
