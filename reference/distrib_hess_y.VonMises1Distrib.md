# von Mises Second Derivative in the Response

Computes \\\partial^2 \ell / \partial y^2\\ in closed form:
\$\$\dfrac{\partial^2 \ell}{\partial y^2} = -\kappa \cos(y - \mu).\$\$
The family is a location family on the circle, so this equals the
pure-direction entry of
[`distrib_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md),
with no sign change. It is negative near the mode and **positive** on
the half of the circle further than a quarter turn from it, the density
having a minimum at the antimode; a periodic density cannot be concave
everywhere.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- y:

  A numeric vector of angles. The expression is evaluated wherever it is
  given, including outside \\\[-\pi, \pi)\\, where the density itself is
  zero.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `kappa` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(kappa))`,
one value per observation.

## See also

[`distrib_grad_y.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.VonMises1Distrib.md)
for the first derivative in the angle,
[`distrib_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md)
for the second derivatives in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, kappa = 2)

# The closed form, written out.
all.equal(distrib_hess_y(d, y, th), -2 * cos(y - 0.5))
#> [1] TRUE

# A location family, so this is the pure-direction entry of the Hessian.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# A central difference of the first derivative reproduces it.
eps <- 1e-5
all.equal((distrib_grad_y(d, y + eps, th) -
           distrib_grad_y(d, y - eps, th)) / (2 * eps),
          distrib_hess_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# Negative at the mode and positive at the antimode: a periodic density is
# not concave everywhere.
distrib_hess_y(d, c(0.5, 0.5 - pi), th)
#> [1] -2  2
```
