# Faa di Bruno Over Set Partitions, on Tables

The partition sum of
[`chain_derivatives`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
taking the inner derivatives as tables rather than fetching them from a
distribution.

## Usage

``` r
chain_assemble(D, inner_params, maps, new_params, order, n)
```

## Arguments

- D:

  A list of length `order`; `D[[k]]` is the inner derivative table of
  order `k`, keyed as
  [`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(inner_params, k)`.

- inner_params:

  The names of the inner coordinates.

- maps:

  Per inner coordinate, a keyed table of the map's partials in the outer
  coordinates, as
  [`reparam_tables`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  returns them.

- new_params:

  The names of the outer parameters.

- order:

  The derivative order, 1 to 4.

- n:

  The length to recycle a constant component to.

## Value

A named list of component vectors.

## Details

Separated so that a family whose own derivatives are easiest to write in
coordinates it does not expose can carry them into the ones it does,
without a second copy of the enumeration.
[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
is the case: its derivatives are compact in the two rates \\(a, c)\\ and
the family is parametrized by \\(\lambda, \alpha)\\, which is a bilinear
map away.

## See also

[`chain_derivatives`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
