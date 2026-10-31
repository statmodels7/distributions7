# Sixth-Order Derivatives by One Tensor Stencil

Computes the sixth-order derivative components, where each component is
a single tensor-product central difference of the highest analytic order
below it, through
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md).

## Usage

``` r
numerical_deriv6(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  accuracy = 4L
)
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  `"parameter"` or `"link"`; the scale on which the differenced quantity
  is read and on which the perturbation is applied.

- accuracy:

  The order of accuracy of the central rules, passed to
  [`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html).
  The default 4 uses four evaluations of the fifth order per parameter
  where it is analytic; 2 uses two.

## Value

A named list of sixth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 6)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

The quantity that is differenced is the fifth order where it is
analytic, which holds for a family that registers its own
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
method and for a wrapper of such a family. Each component is then one
central difference of a fifth-order component in one parameter, and the
whole order costs one stencil per parameter. Where the fifth order is
not analytic, the order that is differenced is the one from which
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md)
differences, so that the sixth order is a stencil of order two (or
higher) on the fourth order, or on the highest order below it that the
family implements. A difference of a difference is never formed.

The step of a difference of order \\r\\ is
\\\varepsilon^{1/(r+a)}\max(1, \|\theta_j\|)\\, with \\a\\ the accuracy,
reduced so that every node stays inside the parameter's domain.

On the link scale the perturbation is applied to \\\eta\\ and the
differenced quantity is read on the link scale, which
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
gives exactly, so the result is the derivative of a link-scale quantity
and needs no sixth derivative of a link.
[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md)
explains why this route is preferred to the composition of the parameter
scale.

## See also

[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md),
the generic;
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
the order below;
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md),
which does the differencing.

## Examples

``` r
numerical_deriv6(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
```
