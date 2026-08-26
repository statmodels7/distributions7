# Gaussian Mixed Second-Response Derivatives

Closed form, and half of it zero. The gaussian's response curvature is
\\\ell^{(yy)} = -1/\sigma^2\\, which carries no location, so
\$\$\frac{\partial^3\ell}{\partial y^2\\\partial\mu} = 0, \qquad
\frac{\partial^3\ell}{\partial y^2\\\partial\sigma} =
\frac{2}{\sigma^3}.\$\$ Neither component varies with the data.

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with `mu` and `sigma`.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors of length `length(y)`, keyed `mu`
and `sigma`, the first exactly zero.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other. \\\ell^{(yy)}\\ is \\\partial^2\ell/\partial y^2\\.

## See also

[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic,
[`distrib_hess_y_hess.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.Gaussian1Distrib.md)
for the next order in \\\theta\\, and
[`distrib_cross2_y.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.StudentT1Distrib.md),
where nothing vanishes.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)
distrib_cross2_y(d, y, theta)
#> $mu
#> [1] 0 0 0
#> 
#> $sigma
#> [1] 0.9103323 0.9103323 0.9103323
#> 

# The one surviving component, written out.
2 / 1.3^3
#> [1] 0.9103323

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::grad(f, c(0.4, 1.3))
#> [1] 0.0000000 0.9103323
```
