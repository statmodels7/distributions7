# Sixth-Order Derivatives

Computes the unique sixth-order partial derivatives of the
log-likelihood with respect to the distribution's parameters, on the
parameter scale or on the link scale.

## Usage

``` r
distrib_deriv6(distrib, y, theta, scale = c("parameter", "link"), ...)
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
  predictors.

- ...:

  Additional arguments passed to the specific method; the default method
  passes `accuracy` to
  [`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md).

## Value

A named list of derivative-component vectors, on the requested scale, in
which each vector has length `length(y)` and the names are those that
[`deriv_names(distrib@params, 6)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives (e.g. `"mu_mu_mu_sigma_sigma_sigma"`).

## Details

The order is needed by a score-driven filter: the derivative of order
\\k\\ of the filtered predictor in the filter's parameters contains the
derivative of order \\k + 1\\ of the log-likelihood in the predictor, so
the fifth derivative of the predictor reads the sixth of the family.

No family computes the sixth order in closed form. The default method,
[`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md),
applies one central difference to the highest order that the family
computes itself: the fifth where
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
is analytic for the family, which holds for every univariate family
shipped with the package that has a closed-form fifth order and for a
wrapper of such a family; otherwise the fourth, or the highest order
below it that the family implements. A difference of a difference is
never formed.

As for
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
there is no `expected` argument.

The route on the link scale depends on where the parameter-scale sixth
order comes from. A family that registers its own method for this
generic (see
[`has_exact_deriv6()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv6.md))
returns the parameter scale, as the methods of orders one to four do,
and the generic composes it with the lower orders by the Faa di Bruno
formula at order six, using the derivatives of the inverse links up to
the sixth. Every term is then exact and the composition is correct up to
rounding. For any other family the sixth order is itself a stencil, and
the stencil is taken in the linear predictors \\\eta\\, on the
link-scale fifth order, which
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
computes exactly from the analytic orders and the derivatives of the
inverse links. Composing a stencil on the parameter scale instead was
measured against 50-digit references and rejected: where a parameter is
small the terms of the Faa di Bruno formula are large and of opposite
sign (for a Poisson mean of 0.05 every link-scale derivative is
\\-\mu\\, while the parameter-scale terms are of order \\120 y/\mu^6\\),
so the error of the stencil is amplified, up to a relative error of 2e-2
against 2e-10 for the stencil in \\\eta\\. The default method reads
`scale` and returns the scale that it is given; a method of a family
does not read it.

## See also

[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the order below,
[`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md)
for the differencing,
[`has_exact_deriv6()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv6.md)
for the predicate that selects the route on the link scale.

## Examples

``` r
distrib_deriv6(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu_mu_mu_mu_mu
#> [1] 0 0 0
#> 
#> $mu_mu_mu_mu_mu_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_mu_mu_sigma_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_mu_sigma_sigma_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_sigma_sigma_sigma_sigma
#> [1] -120 -120 -120
#> 
#> $mu_sigma_sigma_sigma_sigma_sigma
#> [1]  720    0 -720
#> 
#> $sigma_sigma_sigma_sigma_sigma_sigma
#> [1] -2400   120 -2400
#> 
distrib_deriv6(poisson_distrib(), c(0, 2, 5), list(mu = 3), scale = "link")
#> $mu_mu_mu_mu_mu_mu
#> [1] -3 -3 -3
#> 
```
