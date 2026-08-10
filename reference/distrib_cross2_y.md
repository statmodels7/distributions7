# Mixed Second-Response Parameter Derivatives of the Log-Density

Computes \\\partial^3 \ell / \partial y^2\\ \partial \theta_i\\, one
component per parameter, each a vector along `y`.

## Usage

``` r
distrib_cross2_y(distrib, y, theta, scale = c("parameter", "link"), ...)
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

[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
says how curved the log-density is in the response; this says how that
curvature moves with each parameter. It completes the surface a marginal
likelihood needs: where a penalty is \\-\log f(D\beta;\theta)\\, its
Hessian in the coefficients is \\-D'\mathrm{diag}(\ell^{(yy)})D\\ and
the derivative of that in \\\theta\\ is this component placed the same
way.

On the link scale the component for \\\eta_i\\ is the parameter-scale
component multiplied by \\h_i'(\eta_i)\\: the response derivatives are
untouched by a reparametrization of \\\theta\\, so only the first-order
diagonal chain rule enters, exactly as for
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md).

Distributions with a closed form provide it directly; the others fall
back to one central difference of
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
in each parameter (see
[`numerical_cross2_y`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md)).
The reference is the analytic response Hessian, so a family with one
pays for exactly one difference.

## See also

[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)

## Examples

``` r
d <- gaussian1_distrib()
distrib_cross2_y(d, c(-1, 0, 2), list(mu = 0, sigma = 1))
#> $mu
#> [1] 0 0 0
#> 
#> $sigma
#> [1] 2 2 2
#> 
```
