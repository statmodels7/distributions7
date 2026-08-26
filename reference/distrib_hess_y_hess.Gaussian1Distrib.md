# Gaussian Hyperparameter Hessian of the Response Curvature

Closed form, and nearly all of it zero. The gaussian's response
curvature is \\-1/\sigma^2\\, which carries no location at all, so
\$\$\frac{\partial^4\ell}{\partial y^2\partial\sigma^2} =
-\frac{6}{\sigma^4}\$\$ is the only component that does not vanish, and
it does not vary with the data either.

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

A named list of three numeric vectors of length `length(y)`, keyed
`mu_mu`, `sigma_sigma` and `mu_sigma`, of which the first and last are
exactly zero.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other.

## See also

[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the generic,
[`distrib_grad_y_hess.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.Gaussian1Distrib.md)
for the third-order twin, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the first-order column.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)
h <- distrib_hess_y_hess(d, y, theta)
h
#> $mu_mu
#> [1] 0 0 0
#> 
#> $sigma_sigma
#> [1] -2.100767 -2.100767 -2.100767
#> 
#> $mu_sigma
#> [1] 0 0 0
#> 

# One surviving component, constant along y.
c(reported = h$sigma_sigma[1], formula = -6 / 1.3^4)
#>  reported   formula 
#> -2.100767 -2.100767 

# Against a numerical Hessian of the response curvature.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::hessian(f, c(0.4, 1.3))
#>               [,1]          [,2]
#> [1,]  0.000000e+00 -4.174497e-16
#> [2,] -4.174497e-16 -2.100767e+00
```
