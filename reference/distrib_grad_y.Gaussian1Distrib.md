# Gaussian First Derivative in the Response

Computes the first derivative of the Gaussian log-density with respect
to the response, \$\$\dfrac{\partial \ell}{\partial y} = -\dfrac{y -
\mu}{\sigma^2},\$\$ in closed form. The Gaussian is a location family in
\\\mu\\, so the response enters the log-density only through \\y - \mu\\
and this derivative is the negative of the score in \\\mu\\. Quantile
residuals and the mixed derivatives of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
read it.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

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

[`distrib_hess_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gaussian1Distrib.md)
for the second derivative in the response,
[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic and the finite-difference fallback a family without a
closed form takes.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

all.equal(distrib_grad_y(d, y, th), -(y - 0.4) / 1.5^2)
#> [1] TRUE

# A location family: the derivative in the response is minus the score in
# the location.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# It is the derivative of the log-density, so a central difference in y
# reproduces it.
eps <- 1e-5
(distrib_pdf(d, y + eps, th, log = TRUE) -
  distrib_pdf(d, y - eps, th, log = TRUE)) / (2 * eps)
#> [1]  0.71111111  0.04444444 -0.93333333
```
