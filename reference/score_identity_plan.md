# Index Tuples and Partitions of the Score Identity

Returns, for an order, the component names, the parameter index tuple
behind each name and the partitions of the index set into at least two
blocks, which are the ingredients of the score identity of
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## Usage

``` r
score_identity_plan(params, order)
```

## Arguments

- params:

  The parameter names of the family.

- order:

  The derivative order, 2 to 4.

## Value

A list with the elements `nms`, `idx_of` and `parts`.

## Details

The index tuples come from the enumeration that produced the names.
Splitting the names on `"_"` fails for a parameter whose own name
contains one; see
[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md).
At order 2 the names come from
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
which orders the diagonal first, so the tuples come from its inverse and
not from the lexicographic enumeration.

## See also

[`expected_by_bartlett()`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md),
[`score_identity_terms()`](https://statmodels7.github.io/distributions7/reference/score_identity_terms.md)
