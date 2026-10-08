# Bivariate Composition, One Order

The partials of order `order` of \\h(u(x, z))\\ for a scalar outer \\h\\
and a bivariate inner \\u\\, Faa di Bruno written out component by
component: two partials at order one, three at order two, four at order
three and five at order four. The inner partials arrive as a named list
with entries `x`, `z`, `xx`, `xz`, `zz`, `xxx`, `xxz`, `xzz`, `zzz`,
`xxxx`, `xxxz`, `xxzz`, `xzzz`, `zzzz`; missing entries count as zero,
and only the entries of order up to `order` are read. A caller that
needs several orders calls it once per order.

## Usage

``` r
fdb2(h, u, order)
```

## Arguments

- h:

  A list of the outer derivatives `h1` to at least `h_order` at the
  inner value.

- u:

  The named list of inner partials.

- order:

  An integer from 1 to 4.

## Value

A named list of the partials of order `order` of the composition, keyed
as `u` is.

## See also

[`fdb1()`](https://statmodels7.github.io/distributions7/reference/fdb1.md)
for a univariate inner function;
[`gengamma_components()`](https://statmodels7.github.io/distributions7/reference/gengamma_components.md),
[`gpd_components()`](https://statmodels7.github.io/distributions7/reference/gpd_components.md)
and
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md),
which consume it.
