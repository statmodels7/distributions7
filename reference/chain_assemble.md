# Faa di Bruno Over Set Partitions, on Tables

The partition sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
taking the inner derivatives as tables instead of fetching them from a
distribution. With the two split, a family whose derivatives are easiest
in coordinates it does not expose can carry them into the ones it does,
without a second copy of the enumeration.

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

  The length to recycle a constant component to, so that every component
  comes back one value per observation.

## Value

A named list of numeric vectors, keyed as
[`deriv_names(new_params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
is the case the split was made for: its derivatives are compact in the
two rates \\(a, c)\\ and the family is parametrized by \\(\lambda,
\alpha)\\, which is a bilinear map away. The same function also serves
the location-scale cdf derivatives, whose inner coordinate is the
standardized quantile.

The key into an inner table is **built** from the parameter names and
never parsed out of one, so a parameter whose own name contains an
underscore is safe. That is the shape of mistake the package records
having made once, in three separate fallbacks at once.

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
its caller;
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
and
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
the two other consumers;
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
for the enumeration.
