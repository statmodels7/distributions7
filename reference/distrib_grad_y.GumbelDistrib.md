# Gumbel First Derivative in the Response

Computes the first derivative of the Gumbel log-density with respect to
the response, in closed form: \$\$\dfrac{\partial \ell}{\partial y} =
\dfrac{w - 1}{\sigma}, \qquad w = e^{-(y-\mu)/\sigma}.\$\$ The Gumbel is
a location family in \\\mu\\, so the response enters the log-density
only through \\y - \mu\\ and this derivative is the negative of the
score in \\\mu\\. It vanishes at \\y = \mu\\, which is the mode.

Quantile residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

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
one value per observation.

## See also

[`distrib_hess_y.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.GumbelDistrib.md)
for the second derivative in the response,
[`distrib_gradient.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GumbelDistrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)

all.equal(distrib_grad_y(d, y, th), (exp(-y) - 1) / 1)
#> [1] TRUE

# A location family: the derivative in the response is minus the score in
# the location.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Zero at the mode, which is mu itself.
distrib_grad_y(d, 0, th)
#> [1] 0

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-6
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  1.7182818  0.0000000 -0.6321206
```
