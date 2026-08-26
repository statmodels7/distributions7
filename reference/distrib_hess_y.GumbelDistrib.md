# Gumbel Second Derivative in the Response

Computes the second derivative of the Gumbel log-density with respect to
the response, in closed form: \$\$\dfrac{\partial^2 \ell}{\partial y^2}
= -\dfrac{w}{\sigma^2}, \qquad w = e^{-(y-\mu)/\sigma}.\$\$ It is
negative at every observation, \\w\\ being positive, so the log-density
is concave in the response over the whole line. Being a location family,
the Gumbel has the same curvature in the response as in its location,
and this equals the `mu_mu` component of
[`distrib_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md).

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
one value per observation, all strictly negative.

## See also

[`distrib_grad_y.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.GumbelDistrib.md)
for the first derivative in the response,
[`distrib_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)

distrib_hess_y(d, y, th)
#> [1] -2.7182818 -1.0000000 -0.3678794

# A location family: the same curvature as in the location.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# Negative everywhere, so the log-density is concave in the response.
all(distrib_hess_y(d, c(-20, 0, 20), th) < 0)
#> [1] TRUE
```
