# Gaussian Hyperparameter Hessian of the Response Gradient

Closed form. With \\r = y - \mu\\ the response gradient is
\\-r/\sigma^2\\, so differentiating it twice in the parameters gives
\$\$\frac{\partial^3\ell}{\partial y\\\partial\mu^2} = 0, \qquad
\frac{\partial^3\ell}{\partial y\\\partial\mu\\\partial\sigma} =
-\frac{2}{\sigma^3}, \qquad \frac{\partial^3\ell}{\partial
y\\\partial\sigma^2} = -\frac{6r}{\sigma^4}.\$\$ The location pair
vanishes because the response gradient is LINEAR in \\\mu\\, and only
the scale pair varies with the data.

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with `mu` and `sigma`.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors of length `length(y)`, keyed
`mu_mu`, `sigma_sigma` and `mu_sigma`.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other.

## See also

[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
for the generic,
[`distrib_hess_y_hess.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.Gaussian1Distrib.md)
for the fourth-order twin, and
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the first-order column.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)
g <- distrib_grad_y_hess(d, y, theta)
g
#> $mu_mu
#> [1] 0 0 0
#> 
#> $sigma_sigma
#> [1]  2.9410735  0.8403067 -3.3612268
#> 
#> $mu_sigma
#> [1] -0.9103323 -0.9103323 -0.9103323
#> 

# The three formulas written out.
r <- y - 0.4
c(mu_mu = 0, mu_sigma = -2 / 1.3^3, sigma_sigma = -6 * r[1] / 1.3^4)
#>       mu_mu    mu_sigma sigma_sigma 
#>   0.0000000  -0.9103323   2.9410735 

# Against a numerical Hessian of the response gradient.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::hessian(f, c(0.4, 1.3))
#>               [,1]       [,2]
#> [1,] -1.060880e-12 -0.9103323
#> [2,] -9.103323e-01  2.9410735
```
