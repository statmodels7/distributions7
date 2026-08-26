# Weibull Second Derivative in the Response

Computes \\\partial^2 \ell / \partial y^2\\, the second derivative of
the Weibull log-density with respect to the response, in closed form.
With \\u = (y/\mu)^{\sigma}\\, \$\$\dfrac{\partial^2 \ell}{\partial y^2}
= -\dfrac{(\sigma - 1)(1 + \sigma u)}{y^2}.\$\$ The sign is the shape's:
at \\\sigma \> 1\\ the log-density is concave on the whole support, at
\\\sigma \< 1\\ convex, and at \\\sigma = 1\\, the exponential case,
exactly zero, the log-density being linear in \\y\\ there.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of positive observations. At or below zero the result
  is not defined and propagates as `NaN` or `Inf`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma))`,
one value per observation.

## See also

[`distrib_grad_y.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Weibull1Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Weibull1Distrib.md)
for the second derivatives in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)

# The closed form, written out.
u <- (y / 2)^1.5
all.equal(distrib_hess_y(d, y, th), -(1.5 - 1) * (1 + 1.5 * u) / y^2)
#> [1] TRUE

# A central difference of the first derivative reproduces it.
eps <- 1e-5
all.equal((distrib_grad_y(d, y + eps, th) -
           distrib_grad_y(d, y - eps, th)) / (2 * eps),
          distrib_hess_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# At shape one the log-density is linear in y, so the curvature is exactly
# zero; at a smaller shape it turns positive.
rbind(shape_1.0 = distrib_hess_y(d, y, list(mu = 2, sigma = 1)),
      shape_0.6 = distrib_hess_y(d, y, list(mu = 2, sigma = 0.6)))
#>               [,1]      [,2]       [,3]
#> shape_1.0 0.000000 0.0000000 0.00000000
#> shape_0.6 2.017864 0.4004481 0.07845576
```
