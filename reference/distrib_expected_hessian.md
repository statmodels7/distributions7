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

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- method:

  `NULL` (the default) or an object that selects how the expectation is
  obtained. With `NULL`, the closed form of the family is returned where
  the family has one, and
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  (the outer product of the scores, read at each observation without
  taking an expectation) is used otherwise. An object is always used,
  also for a family that has a closed form:
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  (quadrature or summation over the support) or
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
  (Monte Carlo, with its own sample size `nsim`); the last two take
  `identity = "hessian"` for the expectation of the observed Hessian or
  `identity = "score"` for the expectation of minus the outer product of
  the scores. Any other value, a character string included, signals an
  error.
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  compares the methods.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of numeric vectors, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
holding the expected second derivatives, that is minus the Fisher
information.

## Details

On the link scale the first-order term of the chain rule drops out
because the score has zero expectation, so the expected Hessian
transforms as the simple congruence \\\mathrm{diag}(h')\\
\mathbb{E}\[H\]\\ \mathrm{diag}(h')\\ with \\h' = dg^{-1}/d\eta\\.

## See also

[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)

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

# an object is used also by a family with a closed form: the pseudo-Huber
# computes its expected information itself, and the quadrature agrees with it
ph <- pseudohuber_distrib()
th <- list(mu = 0, sigma = 1, nu = 1)
distrib_expected_hessian(ph, 0, th)$mu_mu
#> [1] -0.4545896
distrib_expected_hessian(ph, 0, th, method = expect_integrate())$mu_mu
#> [1] -0.4545896
```
