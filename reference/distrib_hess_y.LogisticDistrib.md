# Logistic Second Derivative in the Response

Computes the second derivative of the logistic log-density with respect
to the response, \$\$\dfrac{\partial^2 \ell}{\partial y^2} =
-\dfrac{1}{2\sigma^2}\\\mathrm{sech}^2\left(\dfrac{z}{2}\right), \qquad
z = \dfrac{y-\mu}{\sigma},\$\$ in closed form. It is negative
everywhere, so the log-density is concave in the response throughout,
and it decays to zero in either tail. Being a location family, the
logistic has the same curvature in the response as in its location, so
this equals the `mu_mu` component of
[`distrib_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md).

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

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

[`distrib_grad_y.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.LogisticDistrib.md)
for the first derivative in the response,
[`distrib_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md)
for the curvature in the parameters,
[`distrib_hess_y.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.CauchyDistrib.md)
for a family that is convex in its tails, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

distrib_hess_y(d, y, th)
#> [1] -0.1693176 -0.2219755 -0.1410532

# A location family: the same curvature as in the location.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# Concave everywhere, unlike the Cauchy, which turns convex in its tails.
all(distrib_hess_y(d, seq(-50, 50, length.out = 1e3), th) < 0)
#> [1] TRUE
any(distrib_hess_y(cauchy_distrib(), seq(-50, 50, length.out = 1e3), th) > 0)
#> [1] TRUE
```
