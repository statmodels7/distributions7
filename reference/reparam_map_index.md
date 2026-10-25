# Index Tables for the Chain Through a Reparametrization

Rewrites the keyed partial tables of a map as lists indexed by an
integer code of the sorted tuple of new coordinates, so that the chain
formulas read a partial without composing a key.

## Usage

``` r
reparam_map_index(maps, p, q, rmax)
```

## Arguments

- maps:

  One keyed table per parent parameter, as
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  returns them.

- p, q:

  The number of parent parameters and of new ones.

- rmax:

  The highest order of partial in `maps`.

## Value

A list with `enc`, the code of a tuple (the empty tuple is 1),
`tup_mat`, the sorted tuples of one length as the rows of a matrix, and
`Jc`, the partial of parent parameter `i` along a tuple of new
coordinates, `0` where the table has no such key.
