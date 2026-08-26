# Inverse Gaussian Second Derivative in the Response

Computes the second derivative of the inverse Gaussian log-density with
respect to the response, in closed form: \$\$\dfrac{\partial^2
\ell}{\partial y^2} = \dfrac{3}{2y^2} - \dfrac{1}{\phi y^3}.\$\$ The
mean drops out, the log-density depending on \\\mu\\ only through terms
linear in \\y\\ and constant in it. The sign changes at \\y =
2/(3\phi)\\: the log-density is concave in the response below that point
and **convex above it**, so an inverse Gaussian log-likelihood is not a
concave function of an observation over the whole support.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of strictly positive observations. The value diverges
  as `y` approaches zero.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. `mu` is not read. `phi` must be
  strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(phi))`, one value per
observation.

## See also

[`distrib_grad_y.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss1Distrib.md)
for the first derivative in the response,
[`distrib_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss1Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)

all.equal(distrib_hess_y(d, y, th), 1.5 / y^2 - 1 / (2 * y^3))
#> [1] TRUE

# The mean does not enter.
identical(distrib_hess_y(d, y, th),
          distrib_hess_y(d, y, list(mu = 300, phi = 2)))
#> [1] TRUE

# Concave below 2/(3 phi) and convex above it.
cut <- 2 / (3 * 2)
c(below = distrib_hess_y(d, cut / 2, th),
  at = distrib_hess_y(d, cut, th),
  above = distrib_hess_y(d, 2 * cut, th))
#>    below       at    above 
#> -54.0000   0.0000   1.6875 
```
