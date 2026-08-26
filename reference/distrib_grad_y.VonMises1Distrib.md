# von Mises First Derivative in the Response

Computes \\\partial \ell / \partial y\\, the derivative of the von Mises
log-density with respect to the angle, in closed form:
\$\$\dfrac{\partial \ell}{\partial y} = -\kappa \sin(y - \mu).\$\$ The
family is a location family on the circle, so this is exactly the
negative of the direction score
[`distrib_gradient.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md)`$mu`.
It is **bounded by** \\\kappa\\ and vanishes at both the mode \\y =
\mu\\ and the antimode \\y = \mu \pm \pi\\, as a periodic density must.

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

[`distrib_hess_y.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.VonMises1Distrib.md)
for the second derivative in the angle,
[`distrib_gradient.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md)
for the derivatives in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, kappa = 2)

# The closed form, written out.
all.equal(distrib_grad_y(d, y, th), -2 * sin(y - 0.5))
#> [1] TRUE

# A location family on the circle, so this is minus the direction score.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# A central difference of the log-density reproduces it.
eps <- 1e-6
all.equal((distrib_pdf(d, y + eps, th, log = TRUE) -
           distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps),
          distrib_grad_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# It vanishes at the mode and at the antimode, and is bounded by kappa.
distrib_grad_y(d, c(0.5, 0.5 - pi), th)
#> [1] 0.000000e+00 2.449294e-16
max(abs(distrib_grad_y(d, seq(-pi, pi, length.out = 401), th)))
#> [1] 1.999993
```
