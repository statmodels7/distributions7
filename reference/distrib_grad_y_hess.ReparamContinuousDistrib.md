# Second-Order Mixed Derivatives of a Reparametrized Distribution

Carries the parent's paired components onto the new parametrization by
the second-order chain rule, through
[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md):
the parent's second-order components multiplied by two first partials of
the map, plus its first-order components multiplied by the map's second
partials. The fourth-order method is the same body reading
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
of the parent instead.

## Arguments

- distrib:

  A `ReparamContinuousDistrib` object, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  A numeric vector of observations.

- theta:

  A named list of the NEW parameters.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per unordered pair of the new
parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

A response derivative does not interact with a reparametrization of
\\\theta\\, so nothing on the \\y\\ side is touched and the whole
correction is the map's. The partials come from
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md),
which a family supplies through `map_derivs` or, failing that, through
one stencil per partial.

## Notation

\\\ell\\ is the log-density of one observation, \\\theta_i\\ a parameter
of the PARENT and \\\alpha_a\\ one of the new parametrization, so that
\\\theta = \theta(\alpha)\\ is the map. \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ are the first and second derivatives in the response.

## See also

[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md),
which does the work,
[`distrib_cross2_y.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.ReparamContinuousDistrib.md)
for the first order in \\\theta\\, and
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
for the wrapper.

## Examples

``` r
# A gaussian obtained in its variance rather than written out.
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma2 = 1.44)

g3 <- distrib_grad_y_hess(d, y, theta)
vapply(g3, function(z) z[1], numeric(1))
#>         mu_mu sigma2_sigma2     mu_sigma2 
#>     0.0000000     0.6697960    -0.4822531 

# Against the family written out by hand, which shares none of the map's
# arithmetic. The map's partials are stencils here, no map_derivs having
# been given, so the two agree to about 1e-09.
max(abs(unlist(g3) -
        unlist(distrib_grad_y_hess(gaussian2_distrib(), y, theta))))
#> [1] 1.844296e-09

# And against a numerical Hessian of the response gradient.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(f, c(0.3, 1.44))
#>               [,1]       [,2]
#> [1,]  9.210739e-13 -0.4822531
#> [2,] -4.822531e-01  0.6697960

# The fourth order takes the same route one derivative up in y.
g4 <- distrib_hess_y_hess(d, y, theta)
max(abs(unlist(g4) -
        unlist(distrib_hess_y_hess(gaussian2_distrib(), y, theta))))
#> [1] 1.676633e-09
```
