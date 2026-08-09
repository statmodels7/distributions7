# Mixed Response-Parameter Derivatives of the Log-Density

Computes the mixed second derivatives \\\partial^2 \ell / \partial y\\
\partial \theta_i\\, one component per parameter, each a vector along
`y`.

## Usage

``` r
distrib_cross_y(distrib, y, theta, scale = c("parameter", "link"), ...)
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

## Details

This is the off-diagonal block of the joint Hessian in \\(y, \theta)\\,
whose diagonal blocks are
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
and
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md).
It is what joint estimation over both arguments needs, and what the
gradient of a profiled objective needs through the implicit function
theorem.

On the link scale the component for \\\eta_i\\ is the parameter-scale
component multiplied by \\h_i'(\eta_i)\\: the response derivative is
untouched by a reparametrization of \\\theta\\, so only the first-order
diagonal chain rule enters, exactly as for
[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md).

Distributions with a closed form provide it directly; the others fall
back to one central difference of
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
in each parameter (see
[`numerical_cross_y`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)).

## See also

[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)

## Examples

``` r
d <- gaussian1_distrib()
distrib_cross_y(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#> $mu
#> [1] 1 1 1
#> 
#> $sigma
#> [1] -2  0  4
#> 
```
