# Analytical Hessian

Computes the analytical observed second derivatives (Hessian matrix) of
the log-likelihood with respect to the distribution's parameters.

## Usage

``` r
distrib_hessian(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- y:

  A numeric vector of observations.

- theta:

  A named list (or named numeric vector) of distribution parameters.
  Each parameter must have length 1 or `length(y)`.

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of numeric vectors, keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`
(e.g. `"mu_sigma"`).

## Examples

``` r
d <- gaussian_distrib()
distrib_hessian(d, c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu
#> [1] -1 -1 -1
#> 
#> $sigma_sigma
#> [1] -2  1 -2
#> 
#> $mu_sigma
#> [1]  2  0 -2
#> 
```
