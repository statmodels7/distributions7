# Elastic-Net Mixed Response-Parameter Derivatives

Computes \\\partial^2\ell/\partial y\\\partial\theta_i\\, one component
per parameter, by differentiating \\\partial\ell/\partial y =
-a\\\mathrm{sgn}(z) - cz\\ in each. With \\z = y-\mu\\ the components
are \\c\\ in the location, \\-\alpha\\\mathrm{sgn}(z) - (1-\alpha)z\\ in
the rate and \\\lambda(z - \mathrm{sgn}(z))\\ in the mixing weight.

The normalizing constant does not appear: it carries no \\y\\, so
differentiating in the response removes it before any parameter is
differentiated.

This is the off-diagonal block of the joint Hessian in \\(y, \theta)\\,
whose diagonals are
[`distrib_hess_y.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.EnetDistrib.md)
and
[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md).
`penalties7` consumes it: a penalty is a negative log-density at the
coefficients, and estimating coefficients and hyperparameters together
needs this block.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations. At `y == mu` the sign is 0 in R, so
  the value returned is one point of the subdifferential.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`.

- scale:

  Either `"parameter"`, the default, or `"link"`. On the link scale the
  chain rule is the gradient's own first-order diagonal one, the
  response derivative not interacting with a reparametrization of
  `theta`; the transformation is applied in the generic's body.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu`, `lambda` and `alpha`, each
of the length of the recycled inputs.

## Notation

\\z = y - \mu\\, \\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\.

## See also

[`distrib_grad_y.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.EnetDistrib.md),
which it differentiates,
[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md)
for the other diagonal block, and
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)
cy <- distrib_cross_y(d, y, th)
names(cy)
#> [1] "mu"     "lambda" "alpha" 

# Against a central difference of the response derivative in a parameter.
eps <- 1e-6
rbind(analytic = cy$lambda,
      numeric = (distrib_grad_y(d, y, list(mu = 0, lambda = 2 + eps,
                                           alpha = 0.5)) -
                 distrib_grad_y(d, y, list(mu = 0, lambda = 2 - eps,
                                           alpha = 0.5))) / (2 * eps))
#>          [,1] [,2] [,3]  [,4]
#> analytic 1.25 0.65 -0.7 -1.55
#> numeric  1.25 0.65 -0.7 -1.55

# The location component is the constant c, and equals minus the second
# response derivative.
c(cross_mu = unique(cy$mu), minus_hess_y = -unique(distrib_hess_y(d, y, th)))
#>     cross_mu minus_hess_y 
#>            1            1 
```
