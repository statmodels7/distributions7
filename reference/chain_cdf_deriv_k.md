# The Chain Rule of Any Order on a Parent's CDF Derivatives

Carries the parent's derivatives of \\F\\ in its own parameters onto the
new ones, for a map given as keyed partial tables.

## Usage

``` r
chain_cdf_deriv_k(parent, q, th_par, maps, new_params, order)
```

## Arguments

- parent:

  The distribution being mapped.

- q:

  A numeric vector of quantiles.

- th_par:

  The parent's parameters at the new ones.

- maps:

  The map's keyed partial tables.

- new_params:

  The new parameter names.

- order:

  The derivative order, 1 to 4.

## Value

A named list of derivative components of \\F\\.

## Details

The parent's tables are taken on the natural scale, `lower.tail = TRUE`
and `log = FALSE`, because the chain rule applies to \\F\\ itself; the
tail and the logarithm are put on afterwards by
[`cdf_scale_k`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md).

## See also

[`chain_assemble`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md),
[`chain_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md)
