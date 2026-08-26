# Lognormal Second-Response and Second-Order Mixed Derivatives

Reads the gaussian's own components at \\t = \log y\\ and divides by a
power of \\y\\. The transformation carries no parameter, so \\t\\ does
not move with \\\theta\\ and only the Jacobian of the RESPONSE
derivatives enters; the parent is `gaussian2`, which carries the same
\\(\mu, \sigma^2)\\, so `theta` passes through unchanged. All three
mixed methods take this route, through
[`lognormal_theta_chain()`](https://statmodels7.github.io/distributions7/reference/lognormal_theta_chain.md).

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- y:

  A numeric vector of observations, strictly positive.

- theta:

  A named list containing `mu` and `sigma2`.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, keyed by parameter for
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
and by parameter pair for
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md).

## Notation

\\t = \log y\\, and \\\ell^{(y)}\\ and \\\ell^{(yy)}\\ are the first and
second derivatives of the log-density in the response.

## See also

[`lognormal_theta_chain()`](https://statmodels7.github.io/distributions7/reference/lognormal_theta_chain.md),
which does the work,
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic, and
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
the general wrapper for a change of variable.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.4, 1.1, 2.3)
theta <- list(mu = 0.2, sigma2 = 0.6)

vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>        mu    sigma2 
#> -10.41667  36.74116 

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::grad(f, c(0.2, 0.6))
#> [1] -10.41667  36.74116

# The second order in theta, and its numerical Hessian.
g3 <- distrib_grad_y_hess(d, y, theta)
vapply(g3, function(z) z[1], numeric(1))
#>         mu_mu sigma2_sigma2     mu_sigma2 
#>      0.000000     25.840063     -6.944444 
h <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(h, c(0.2, 0.6))
#>               [,1]      [,2]
#> [1,] -1.800539e-10 -6.944444
#> [2,] -6.944444e+00 25.840063
```
