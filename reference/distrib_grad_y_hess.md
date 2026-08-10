# Hyperparameter Hessians of the Response Derivatives

`distrib_grad_y_hess()` computes \\\partial^3 \ell / \partial y\\
\partial\theta_i \partial\theta_j\\ and `distrib_hess_y_hess()` computes
\\\partial^4 \ell / \partial y^2\\ \partial\theta_i \partial\theta_j\\,
one component per unordered pair of parameters, each a vector along `y`.

## Usage

``` r
distrib_grad_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)

distrib_hess_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- y:

  A numeric vector of observations.

- theta:

  A named list (or named numeric vector) of distribution parameters.

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- ...:

  Additional arguments passed to the specific method.

## Value

A named list with one numeric vector per parameter pair, keyed by
`hess_names(distrib@params)`.

## Details

These are the second-order column of the mixed grid whose first-order
one is
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
and
[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md):
how the response gradient and the response curvature CURVE in the
parameters. A marginal likelihood needs them to differentiate its
Laplace approximation twice, the penalty being a negative log-density
evaluated at the coefficients.

The components are keyed by
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
the same enumeration
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
uses, so a consumer that looks a pair up finds it under the name it
already knows.

On the link scale each component is multiplied by \\h_i'(\eta_i)
h_j'(\eta_j)\\ and, on a diagonal pair, gains the second-order term
\\h_i''(\eta_i)\\ times the corresponding first-order component: the
response derivatives are untouched by a reparametrization of \\\theta\\,
so this is the ordinary diagonal chain rule at second order, exactly as
for
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md).

Distributions with closed forms provide them; the others take one
central difference of the analytic first-order quantity in each
parameter (see
[`numerical_theta2_y`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)).
The two differences act on different parameters off the diagonal, so
they compose into a single mixed stencil rather than the nested
differencing of one variable the package forbids.

## See also

[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)

## Examples

``` r
d <- gaussian1_distrib()
distrib_grad_y_hess(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#> $mu_mu
#> [1] 0 0 0
#> 
#> $sigma_sigma
#> [1]   6   0 -12
#> 
#> $mu_sigma
#> [1] -2 -2 -2
#> 
distrib_hess_y_hess(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#> $mu_mu
#> [1] 0 0 0
#> 
#> $sigma_sigma
#> [1] -6 -6 -6
#> 
#> $mu_sigma
#> [1] 0 0 0
#> 
```
