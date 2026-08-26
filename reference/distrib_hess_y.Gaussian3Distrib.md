# Gaussian Second Derivative in the Response, Mean and Precision

Computes the second derivative of the Gaussian log-density with respect
to the response, \$\$\dfrac{\partial^2 \ell}{\partial y^2} = -\tau,\$\$
in closed form. It does not depend on \\y\\ or on \\\mu\\, so the value
is constant within a parameter setting and is recycled to the length of
`y`. Being a location family, the Gaussian has the same curvature in the
response as in its location, and this equals the `mu_mu` component of
[`distrib_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian3Distrib.md).
In this parametrization that curvature is the parameter itself, which is
the sense in which the precision measures how sharply the log-density
peaks.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `y`. `mu` is not read. `tau` must be
  strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `length(y)`, every entry \\-\tau\\.

## See also

[`distrib_grad_y.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gaussian3Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian3Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, tau = 0.25)

distrib_hess_y(d, y, th)
#> [1] -0.25 -0.25 -0.25

# A location family: the same curvature as in the location, and here that
# curvature is minus the parameter.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# Negative everywhere, so the log-density is concave in the response.
all(distrib_hess_y(d, y, th) < 0)
#> [1] TRUE
```
