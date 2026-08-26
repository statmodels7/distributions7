# Cauchy Second Derivative in the Response

Computes the second derivative of the Cauchy log-density with respect to
the response, \$\$\dfrac{\partial^2 \ell}{\partial y^2} = \dfrac{2(r^2 -
\sigma^2)}{(\sigma^2 + r^2)^2}, \qquad r = y - \mu,\$\$ in closed form.
Being a location family, the Cauchy has the same curvature in the
response as in its location, so this equals the `mu_mu` component of
[`distrib_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md).
It is negative only for \\\|r\| \< \sigma\\: the log-density is concave
near the location and convex in the tails, which is how heavy tails
appear on the log scale.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

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

[`distrib_grad_y.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.CauchyDistrib.md)
for the first derivative in the response,
[`distrib_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

distrib_hess_y(d, y, th)
#> [1]  0.02679795 -0.87712429  0.09739469

# A location family: the same curvature as in the location.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# Concave within one scale unit of the location, convex beyond it.
z <- c(0, 1, 1.5, 2, 5)
data.frame(r = z * 1.5,
           hess_y = distrib_hess_y(d, 0.4 + z * 1.5, th),
           concave = distrib_hess_y(d, 0.4 + z * 1.5, th) < 0)
#>      r      hess_y concave
#> 1 0.00 -0.88888889    TRUE
#> 2 1.50  0.00000000   FALSE
#> 3 2.25  0.10519395   FALSE
#> 4 3.00  0.10666667   FALSE
#> 5 7.50  0.03155819   FALSE
```
