# Analytical Gradient

Computes the analytical first derivatives of the log-likelihood with
respect to the distribution's parameters.

## Usage

``` r
distrib_gradient(distrib, y, theta, scale = c("parameter", "link"), ...)
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

A named list with one numeric vector per parameter, keyed by
`distrib@params`.

## Examples

``` r
d <- gaussian1_distrib()
distrib_gradient(d, c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu
#> [1] -1  0  1
#> 
#> $sigma
#> [1]  0 -1  0
#> 

# the same score with respect to the unconstrained parameters
distrib_gradient(d, c(-1, 0, 1), list(mu = 0, sigma = 1), scale = "link")
#> $mu
#> [1] -1  0  1
#> 
#> $sigma
#> [1]  0 -1  0
#> 
```
