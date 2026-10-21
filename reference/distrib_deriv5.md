# Fifth-Order Derivatives

Computes the unique fifth-order partial derivatives of the
log-likelihood with respect to the distribution's parameters, on the
parameter scale or on the link scale.

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
  predictors. In the numerical default the link scale is obtained by
  differentiating the order-4 component already on the link scale, which
  needs neither the fifth derivative of a link nor an order-5 entry in
  [`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md).
  A family with a closed-form fifth order composes it with the lower
  orders through
  [`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md),
  which uses both.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list of derivative-component vectors, in which each vector has
length `length(y)` and the names are those that
[`deriv_names(distrib@params, 5)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives (e.g. `"mu_mu_sigma_sigma_sigma"`).

## Details

The order is needed because each order of differentiating the predictor
of a score-driven filter through its own recursion requires one more
order of the family: the curvature reaches the third order, the
directional third derivative the fourth, and the outer Hessian of a
model with such a term the fifth. The fifth order is therefore provided
by the family, like the lower orders.

The families that register their own method
([`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md),
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md),
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md),
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md),
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md),
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md),
[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md),
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md),
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md),
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md),
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md),
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md),
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md),
[`pseudohuber2_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber2_distrib.md),
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md),
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md),
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
and
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md))
compute the fifth order in closed form with compiled kernels. Every
other family reaches
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
which applies one central difference to the fourth order. Where the
fourth order is analytic this is an ordinary numerical derivative; where
the fourth order is itself a fallback the result is a difference of a
difference, with the accuracy described on
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md).

Unlike the generics of lower order this one carries no `expected`
argument. The criterion of a model with a filter reads the observed
information (a filter has no expected one), and an expectation of the
fifth order would sit on top of a numerical quantity.
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
is generic in the order, so the expected fifth order can be added when a
caller needs it.

## See also

[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order below,
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md)
for the differencing, and
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
for the predicate that reports whether the quantity that is differenced
is itself analytic.

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
