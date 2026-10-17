# Faa di Bruno Over Set Partitions, on Tables

Computes the partition sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
taking the inner derivatives as tables instead of fetching them from a
distribution. Because the two functions are separate, a family whose
derivatives are easiest in coordinates that it does not expose can carry
them into the coordinates that it does expose, without a second copy of
the enumeration.

## Usage

``` r
chain_assemble(D, inner_params, maps, new_params, order, n)
```

## Arguments

- D:

  A list of length `order`. `D[[k]]` is the inner derivative table of
  order \\k\\, keyed as
  [`deriv_names(inner_params, k)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).
  `D[[1]]` may be `NULL` where the first-order term is known to vanish.

- inner_params:

  A character vector naming the inner coordinates.

- maps:

  Per inner coordinate, a keyed table of the map's partials in the outer
  coordinates, as
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  returns them.

- new_params:

  A character vector naming the outer parameters.

- order:

  The derivative order, 1 to 4.

- n:

  The length to which a constant component is recycled, so that every
  component comes back one value per observation.

## Value

A named list of numeric vectors, keyed as
[`deriv_names(new_params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
is the case for which the split was made: its derivatives are compact in
the two rates \\(a, c)\\, and the family is parametrized by \\(\lambda,
\alpha)\\, which are related to \\(a, c)\\ by a bilinear map. The same
function also serves the location-scale cdf derivatives, whose inner
coordinate is the standardized quantile.

The key into an inner table is built from the parameter names and never
parsed out of one, so a parameter whose own name contains an underscore
is safe.

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
its caller;
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md),
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
and the mapped cdf derivatives of R/cdf_mapped_higher.R, further
consumers;
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
for the enumeration.
