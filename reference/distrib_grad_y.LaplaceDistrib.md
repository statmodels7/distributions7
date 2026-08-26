# Laplace First Derivative in the Response

Computes the first derivative of the Laplace log-density with respect to
the response, \$\$\dfrac{\partial \ell}{\partial y} =
-\dfrac{\mathrm{sign}(y - \mu)}{\sigma},\$\$ in closed form. The Laplace
is a location family in \\\mu\\, so this is the negative of the score in
\\\mu\\. It takes only three values, and at \\y = \mu\\ exactly the
method returns 0, the midpoint of the subdifferential; the derivative
does not exist there, the log-density having a corner.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma))`,
taking the values \\-1/\sigma\\, 0 and \\1/\sigma\\ only.

## See also

[`distrib_hess_y.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LaplaceDistrib.md)
for the second derivative, which is zero;
[`distrib_gradient.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md)
for the score in the parameters; and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

all.equal(distrib_grad_y(d, y, th), -sign(y - 0.4) / 1.5)
#> [1] TRUE

# A location family: minus the score in the location.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Three values only, and 0 at the kink itself.
distrib_grad_y(d, 0.4 + c(-100, -1e-9, 0, 1e-9, 100), th)
#> [1]  0.6666667  0.6666667  0.0000000 -0.6666667 -0.6666667
```
