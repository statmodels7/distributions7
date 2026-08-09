# Analytical Expected Hessian

Computes the analytical expected second derivatives of the
log-likelihood with respect to the distribution's parameters.

## Usage

``` r
distrib_expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000,
  ...
)
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

- approx:

  How the expectation is approximated **when the distribution has no
  closed-form expected Hessian**; ignored otherwise. One of `"bartlett"`
  (default, equivalently `"opg"`: the outer product of the score),
  `"integrate"` (quadrature/summation of the observed Hessian) or `"mc"`
  (Monte Carlo). See
  [`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  for the accuracy/speed trade-offs.

- nsim:

  Monte Carlo sample size used when `approx = "mc"`. Defaults to 10000.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of numeric vectors, keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`,
holding the expected second derivatives, that is minus the Fisher
information.

## Details

On the link scale the first-order term of the chain rule drops out
because the score has zero expectation, so the expected Hessian
transforms as the simple congruence \\\mathrm{diag}(h')\\
\mathbb{E}\[H\]\\ \mathrm{diag}(h')\\ with \\h' = dg^{-1}/d\eta\\.

## See also

[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)

## Examples

``` r
d <- gaussian1_distrib()
distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1))
#> $mu_mu
#> [1] -1
#> 
#> $sigma_sigma
#> [1] -2
#> 
#> $mu_sigma
#> [1] 0
#> 

# a family with no closed form uses the strategy named by 'approx'
distrib_expected_hessian(
  pseudohuber_distrib(), 0, list(mu = 0, sigma = 1, nu = 1),
  approx = "integrate"
)
#> $mu_mu
#> [1] -0.4545896
#> 
#> $sigma_sigma
#> [1] -1.244894
#> 
#> $nu_nu
#> [1] -0.01403316
#> 
#> $mu_sigma
#> [1] 0
#> 
#> $mu_nu
#> [1] 0
#> 
#> $sigma_nu
#> [1] -0.1224472
#> 
```
