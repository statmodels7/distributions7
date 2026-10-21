# CDF Derivatives of a Family Written as a Map of Another

This is the
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
body that every mapped family shares. The body tests whether the
parent's cdf derivatives are exact at every order up to the one wanted.
If they are, it calls
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md);
if they are not, it calls
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
on the new family's own cdf.

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

  The distribution onto which the mapped family maps.

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

  A single logical; `TRUE` selects the lower tail and `FALSE` the upper
  tail.

- log:

  A single logical; `TRUE` selects the derivatives of the log
  probability and `FALSE` those of the probability.

## Value

A named list of numeric vectors on the requested tail and scale: one per
new parameter at order 1, and one per
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
component at order 2.

## The exactness test

[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
is applied to the parent at order 1 and, for a Hessian, at order 2 as
well. The chain is exact only if what it carries is exact; if it is
taken over a differenced parent, it adds rounding to rounding and costs
more than differencing the child's cdf once. The test is made at each
call, so a parent that gains a closed order moves its mapped families
onto the chain without any edit here.

## See also

[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
for the exact route;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the fallback;
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
for the gate.
