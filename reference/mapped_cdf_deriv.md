# CDF Derivatives Through a Map, With the Fallback as a Guard

The body shared by the families that obtain their cdf derivatives from a
parent's: the chain rule of
[`chain_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
when the parent has a closed form at that order, and the
finite-difference fallback otherwise.

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

  The distribution in the new parametrization.

- parent:

  The parent distribution.

- th_par:

  The parent's parameters at the new ones.

- maps:

  The map's keyed partial tables.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the new parameters.

- order:

  The derivative order, 1 or 2.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

## Value

A named list of derivative component vectors.

## See also

[`chain_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
