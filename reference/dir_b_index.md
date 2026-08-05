# The Position of a Second-Derivative Block

Locates the \\(k, l)\\ entry in the list of second derivatives a
parameters7 parameter returns, whose elements are keyed by unordered
tuple rather than by position.

## Usage

``` r
dir_b_index(idx, k, l)
```

## Arguments

- idx:

  The tuple index list, as returned by
  [`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html).

- k, l:

  The two free-value positions.

## Value

An integer position into the second-derivative list.

## See also

[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
