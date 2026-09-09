# Fifth-Order Derivatives

The unique fifth-order partial derivatives of the log-likelihood with
respect to the distribution's parameters, on the parameter scale or on
the link scale.

## Usage

``` r
distrib_deriv5(distrib, y, theta, scale = c("parameter", "link"), ...)
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

  `"parameter"` (the default) for derivatives with respect to the
  parameters, `"link"` for derivatives with respect to the unconstrained
  predictors. The link scale is obtained by differentiating the order-4
  component **already on the link scale**, so it needs neither the fifth
  derivative of a link nor an order-5 entry in
  [`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md).

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of derivative-component vectors, each of length
`length(y)`, keyed as
[`deriv_names(distrib@params, 5)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them (e.g. `"mu_mu_sigma_sigma_sigma"`).

## Details

The order exists because each order of differentiating a score-driven
filter's predictor through its own recursion draws in one more order of
the family: the curvature reaches the third, the directional third
derivative the fourth, and the outer Hessian of a model carrying such a
term the fifth. Writing it inside the term that needs it would be the
private helper this package exists to abolish, so it is a family
quantity like any other.

No family computes it in closed form yet. Every one reaches
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
which applies **one** central difference to the analytic fourth order –
not a difference of a difference, which is the rule the whole derivative
surface obeys. Against Richardson extrapolation on the same analytic
fourth that is worth about 1e-10, so the numerical fifth is accurate
enough to build on rather than a placeholder for the architecture.

Unlike its siblings this generic carries no `expected` argument. The
criterion of a model with a filter reads the observed information – a
filter has no expected one – and an expectation of the fifth order would
sit on top of a numerical quantity, which is two approximations deep.
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
is generic in the order, so the expected fifth can be added the day
something asks for it.

## See also

[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order below,
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md)
for the differencing, and
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
for the predicate that says whether the quantity being differenced is
itself analytic.

## Examples

``` r
distrib_deriv5(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu_mu_mu_mu
#> [1] 0 0 0
#> 
#> $mu_mu_mu_mu_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_mu_sigma_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_sigma_sigma_sigma
#> [1] 24 24 24
#> 
#> $mu_sigma_sigma_sigma_sigma
#> [1] -120    0  120
#> 
#> $sigma_sigma_sigma_sigma_sigma
#> [1] 336 -24 336
#> 
```
