# Hyperparameter Hessian of the Response Gradient

Computes \\\partial^3 \ell / \partial y\\ \partial\theta_i
\partial\theta_j\\, one component per unordered pair of parameters, each
a vector along `y`. It says how the response GRADIENT curves in the
parameters, and it is the third-order entry of the mixed grid below.

## Usage

``` r
distrib_grad_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A distribution object inheriting from `continuous_distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list, or named numeric vector, of distribution parameters.
  Aligned by the generic before dispatch.

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- ...:

  Passed to the method.

## Value

A named list with one numeric vector per unordered pair of parameters,
each as long as `y`, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

These are the second-order column of the mixed grid whose first-order
column is
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md):

|  |  |  |  |
|----|----|----|----|
|  | \\\partial^0/\partial\theta^0\\ | \\\partial^1/\partial\theta^1\\ | \\\partial^2/\partial\theta^2\\ |
| \\\partial/\partial y\\ | [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md) | [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md) | `distrib_grad_y_hess()` |
| \\\partial^2/\partial y^2\\ | [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md) | [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md) | [`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md) |

The right-hand column is the one a marginal criterion's SECOND
derivative asks of a penalty. A penalty is a negative log-density read
at the coefficients, so \\\partial^2\rho/\partial\beta^2\\ carries the
density's response curvature and the two mixed derivatives above carry
\\\partial^3\rho/\partial\beta\\\partial\theta^2\\ and
\\\partial^4\rho/\partial\beta^2\partial\theta^2\\.

The components are keyed by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
the same enumeration
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
uses, so a consumer looking a pair up finds it under the name it already
knows.

## The link scale

Each component is multiplied by \\h_i'(\eta_i) h_j'(\eta_j)\\ and, on a
diagonal pair, gains \\h_i''(\eta_i)\\ times the corresponding
first-order component from
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md).
The response derivative is untouched by a reparametrization of
\\\theta\\, so this is the ordinary diagonal chain rule at second order,
exactly as for
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md).

## Where the numbers come from

A distribution with a closed form provides it. The rest take one central
difference of the analytic
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
in each parameter, through
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md).
Off the diagonal the two differences act on DIFFERENT parameters, so
they compose into one mixed stencil rather than the nested differencing
of a single variable the package forbids.

Continuous distributions only, as with every other response derivative:
a discrete family has no method and the call raises with its class
named.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other.

## See also

[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the fourth-order twin,
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the first-order column, and
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)
for the fallback.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)
distrib_grad_y_hess(d, y, theta)
#> $mu_mu
#> [1] 0 0 0
#> 
#> $sigma_sigma
#> [1]  2.9410735  0.8403067 -3.3612268
#> 
#> $mu_sigma
#> [1] -0.9103323 -0.9103323 -0.9103323
#> 

# Closed form: the response gradient is -r / sigma^2, so the location pair
# vanishes and the other two are elementary.
r <- y - 0.4
c(mu_sigma = -2 / 1.3^3, sigma_sigma = -6 * r[1] / 1.3^4)
#>    mu_sigma sigma_sigma 
#>  -0.9103323   2.9410735 

# Against a numerical Hessian of the response gradient in the parameters.
g <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::hessian(g, c(0.4, 1.3))
#>               [,1]       [,2]
#> [1,] -1.060880e-12 -0.9103323
#> [2,] -9.103323e-01  2.9410735

# On the link scale, with sigma on a log link: h' = h'' = sigma.
cy <- distrib_cross_y(d, y, theta)
distrib_grad_y_hess(d, y, theta, scale = "link")$sigma_sigma
#> [1]  3.3136095  0.9467456 -3.7869822
distrib_grad_y_hess(d, y, theta)$sigma_sigma * 1.3^2 + cy$sigma * 1.3
#> [1]  3.3136095  0.9467456 -3.7869822

# A discrete family has no method at all.
try(distrib_grad_y_hess(poisson_distrib(), 1:3, list(mu = 2)))
#> Error : Can't find method for `distrib_grad_y_hess(<distributions7::PoissonDistrib>)`.
```
