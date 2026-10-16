# Fifth-Order Derivatives by One Tensor Stencil

The fifth-order derivative components, each a single tensor-product
central difference of the highest analytic order below it, through
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md).

## Usage

``` r
numerical_deriv5(
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
  The default 4 uses four evaluations of the fourth order per parameter
  where it is analytic; 2 uses two.

## Value

A named list of fifth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 5)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

The quantity differenced is the fourth order where it is analytic, as
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
reads it, and each component is then one central difference of a
fourth-order component in one parameter: the component of a
non-decreasing five-tuple is the four-tuple left when one index is
dropped, and the whole order costs one stencil per parameter. Where the
fourth order is not analytic, a family that registers no wrapper reads
the highest order it implements itself, as
[`analytic_order()`](https://statmodels7.github.io/distributions7/reference/analytic_order.md)
gives it, and a wrapper of a parent without an analytic fourth order
differences its log-density, because its lower orders are partition sums
over the parent's fallbacks. A difference of a difference is never
formed.

The step of a difference of order \\r\\ is
\\\varepsilon^{1/(r+a)}\max(1, \|\theta_j\|)\\, with \\a\\ the accuracy,
reduced so that every node stays inside the parameter's domain. On the
link scale the perturbation is applied to \\\eta\\ and the differenced
quantity is read on the link scale as well, so the result is the
derivative of a link-scale quantity. This route requires neither the
fifth derivative of a link nor an order-5 entry in
[`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md),
both of which the chain rule route uses.

## See also

[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
the generic;
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
the order below;
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md),
which does the differencing.

## Examples

``` r
numerical_deriv5(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
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
