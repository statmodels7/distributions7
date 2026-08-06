# The Chain Rule on a Parent's CDF Derivatives

Carries the derivatives of the parent's distribution function onto a new
parametrization, given the map's partial derivatives.

## Usage

``` r
chain_cdf_deriv(parent, q, th_par, maps, new_params, order)
```

## Arguments

- parent:

  The distribution being reparametrized.

- q:

  A numeric vector of quantiles.

- th_par:

  The parent's parameters at the new ones.

- maps:

  The map's keyed partial tables.

- new_params:

  The new parameter names.

- order:

  The derivative order, 1 or 2.

## Value

A named list of derivative components of \\F\\.

## Details

The map derivatives arrive in the keyed form
[`reparam_map_derivs`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
produces, so a missing key is an exact zero, and the result is exact
whenever the parent's own cdf derivatives are. It is the caller's
business to check that they are: applied to a parent whose derivatives
are themselves differenced, the chain adds an exact transformation to an
approximate quantity and buys nothing over differencing the new cdf
directly.

## See also

[`cdf_tail_scale`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md),
[`chain_derivatives`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
