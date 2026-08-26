# Second-Response Mixed Derivatives of a Reparametrized Distribution

Carries the parent's block onto the new parametrization by the
FIRST-order chain rule, exactly as
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
is carried: a derivative in the response does not interact with a
reparametrization of the parameters, so only the map's first partials
enter and its second ones never appear.

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

A named list with one numeric vector per new parameter, keyed by
`distrib@params`.

## Notation

\\\ell\\ is the log-density of one observation, \\\theta_i\\ a parameter
of the PARENT and \\\alpha_a\\ one of the new parametrization, so that
\\\theta = \theta(\alpha)\\ is the map. \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ are the first and second derivatives in the response.

## See also

[`mapped_cross2_y()`](https://statmodels7.github.io/distributions7/reference/mapped_cross2_y.md),
which does the work,
[`distrib_grad_y_hess.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.ReparamContinuousDistrib.md)
for the second order in \\\theta\\, and
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
for the wrapper.

## Examples

``` r
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

vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>        mu    sigma2 
#> 0.0000000 0.4822531 

# Against the family written out by hand, and against a numerical
# derivative of the response Hessian.
max(abs(unlist(distrib_cross2_y(d, y, theta)) -
        unlist(distrib_cross2_y(gaussian2_distrib(), y, theta))))
#> [1] 6.324885e-12
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::grad(f, c(0.3, 1.44))
#> [1] 0.0000000 0.4822531
```
