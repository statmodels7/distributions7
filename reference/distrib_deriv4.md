# Analytical Fourth-Order Derivatives

Computes the unique fourth-order partial derivatives of the
log-likelihood with respect to the distribution's parameters.
Distributions with a closed-form implementation provide it directly (in
C++); the others fall back to finite differences of the Hessian (see
[`numerical_deriv4`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)).

## Usage

``` r
distrib_deriv4(
  distrib,
  y,
  theta,
  expected = FALSE,
  scale = c("parameter", "link"),
  approx = c("integrate", "bartlett", "mc", "opg"),
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

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives
  (\\\mathbb{E}\[\partial^4 \ell\]\\) instead of the observed ones.
  Defaults to `FALSE`.

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- approx:

  How the expectation is approximated when `expected = TRUE` and the
  distribution has no closed-form expected fourth derivatives; ignored
  otherwise. One of `"integrate"` (default), `"bartlett"` (equivalently
  `"opg"`) or `"mc"`. See
  [`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- nsim:

  Monte Carlo sample size used when `approx = "mc"`. Defaults to 10000.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of derivative-component vectors, keyed as in
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(distrib@params, 4)`
(e.g. `"mu_mu_sigma_sigma"`).

## Examples

``` r
distrib_deriv4(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu_mu_mu
#> [1] 0 0 0
#> 
#> $mu_mu_mu_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_sigma_sigma
#> [1] -6 -6 -6
#> 
#> $mu_sigma_sigma_sigma
#> [1]  24   0 -24
#> 
#> $sigma_sigma_sigma_sigma
#> [1] -54   6 -54
#> 
```
