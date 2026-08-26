# Laplace First Derivative in the Response, Rate Parametrization

Computes the first derivative of the Laplace log-density with respect to
the response, \$\$\dfrac{\partial \ell}{\partial y} =
-\lambda\\\mathrm{sign}(y - \mu),\$\$ in closed form. The family is a
location family in \\\mu\\, so this is the negative of the score in
\\\mu\\. It takes only the three values \\-\lambda\\, 0 and \\\lambda\\;
at \\y = \mu\\ exactly, `sign(0)` is 0 and the method returns 0, the
midpoint of the subdifferential, the derivative not existing there.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `lambda` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(lambda))`,
taking the values \\-\lambda\\, 0 and \\\lambda\\ only.

## See also

[`distrib_hess_y.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Laplace2Distrib.md)
for the second derivative, which is zero;
[`distrib_grad_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LaplaceDistrib.md)
for the scale parametrization;
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, lambda = 2)

all.equal(distrib_grad_y(d, y, th), -2 * sign(y - 0.4))
#> [1] TRUE

# A location family: minus the score in the location.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Three values only, and 0 at the kink itself.
distrib_grad_y(d, 0.4 + c(-100, -1e-9, 0, 1e-9, 100), th)
#> [1]  2  2  0 -2 -2
```
