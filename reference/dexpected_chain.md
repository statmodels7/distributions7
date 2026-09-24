# The Chain Rule for the Expected Information's Derivatives

The arithmetic of
[`reparam_dexpected()`](https://statmodels7.github.io/distributions7/reference/reparam_dexpected.md),
given the parent's quantities rather than the parent: \\E\\, \\\partial
E\\ and, at order 2, \\\partial^2 E\\ on the parent's scale, and the
map's partials as keyed tables. It serves a reparametrized family and
any family written as a map of another's coordinates, such as
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
over the shapes.

## Usage

``` r
dexpected_chain(Pp, Pn, E, d1, d2, maps, order, n)
```

## Arguments

- Pp, Pn:

  The parent's parameter names and the new ones.

- E, d1, d2:

  Named lists keyed as
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
  and
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
  on the parent's parameters; `d2` is read only at order 2.

- maps:

  One keyed table per parent parameter, as
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  returns them: the key is the sorted, comma-joined indices of the new
  parameters differentiated, and a missing key is an exact zero.

- order:

  `1L` or `2L`.

- n:

  The number of observations, to which every component is recycled.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
on the new parameters.

## See also

[`reparam_dexpected()`](https://statmodels7.github.io/distributions7/reference/reparam_dexpected.md)
