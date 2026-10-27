# Analytical Third-Order Derivatives

Computes the unique third-order partial derivatives of the
log-likelihood with respect to the distribution's parameters.
Distributions with a closed-form implementation provide it directly (in
C++); the others fall back to finite differences of the highest order
that the family implements itself (see
[`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md)).

## Usage

``` r
distrib_deriv3(
  distrib,
  y,
  theta,
  expected = FALSE,
  scale = c("parameter", "link"),
  method = NULL,
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

  Logical; if `TRUE`, returns the expected third derivatives
  (\\\mathbb{E}\[\partial^3 \ell\]\\) instead of the observed ones.
  Defaults to `FALSE`.

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- method:

  `NULL` (the default) or an object that selects how the expectation is
  obtained; it is read only when `expected = TRUE`. With `NULL`, the
  closed form of the family is returned where the family has one, and
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  (with `identity = "hessian"`) is used otherwise. An object is always
  used, also for a family that has a closed form:
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  (the empirical sum of products of the score identity, with no
  expectation),
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  or
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md),
  the last two with `identity = "hessian"` or `identity = "score"`. Any
  other value, a character string included, signals an error.
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  compares the methods.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of derivative-component vectors, keyed as in
[`deriv_names(distrib@params, 3)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
(e.g. `"mu_mu_sigma"`).

## See also

[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)

## Examples

``` r
distrib_deriv3(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu_mu
#> [1] 0 0 0
#> 
#> $mu_mu_sigma
#> [1] 2 2 2
#> 
#> $mu_sigma_sigma
#> [1] -6  0  6
#> 
#> $sigma_sigma_sigma
#> [1] 10 -2 10
#> 
```
