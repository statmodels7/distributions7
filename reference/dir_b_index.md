# The Position of a Second-Derivative Block

Returns the position of the \\(k, l)\\ entry in the list of second
derivatives a `parameters7` parameter supplies. That list is keyed by
unordered tuple rather than by position, so a caller holding two
free-value indices has to search for the entry naming them, in either
order.

## Usage

``` r
dir_b_index(idx, k, l)
```

## Arguments

- idx:

  The tuple index list, as returned by
  [`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html)
  at order 2.

- k, l:

  The two free-value positions, in either order.

## Value

A single integer, the position into the second-derivative list. `NA` if
no tuple matches, which cannot happen for indices inside the parameter's
own range.

## See also

[`dir_parts()`](https://statmodels7.github.io/distributions7/reference/dir_parts.md),
which supplies both the list and the tuples, and
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
for the family.
