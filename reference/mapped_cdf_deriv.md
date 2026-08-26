# CDF Derivatives of a Family Written as a Map of Another

The
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
body every mapped family shares. It asks whether the parent's cdf
derivatives are exact at every order up to the one wanted, takes
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
if they are, and falls back to
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
on the new family's own cdf if they are not.

## Usage

``` r
mapped_cdf_deriv(
  distrib,
  parent,
  th_par,
  maps,
  q,
  theta,
  order,
  lower.tail,
  log
)
```

## Arguments

- distrib:

  The mapped family, whose cdf and parameter names are read.

- parent:

  The distribution it maps onto.

- th_par:

  The parent's parameters, evaluated at the new ones.

- maps:

  The map's keyed partial tables.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters.

- order:

  The derivative order, 1 or 2.

- lower.tail:

  Is the lower tail wanted? A single logical.

- log:

  Are derivatives of the log probability wanted? A single logical.

## Value

A named list of numeric vectors on the requested tail and scale: one per
new parameter at order 1, and one per
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
component at order 2.

## The gate, and why it is a gate

[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
is asked of the parent at order 1 and, for a Hessian, at order 2 as
well. The chain is exact only if what it carries is exact; taken over a
differenced parent it would add rounding to rounding and cost more than
differencing the child's cdf once. The gate is also what keeps a
family's page honest as the parent improves: when the inverse Gaussian
gained a closed second order, its second parametrization stopped
differencing without any edit here.

## See also

[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the exact route;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the fallback;
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
for the gate.
