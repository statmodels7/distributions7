# The Second Derivative of the Covariance, by Position

`param_d2` keyed by the pair of free values rather than by the string
the structure names it with, the key being CONSTRUCTED from the sorted
pair and never parsed out of a name.

## Usage

``` r
mvg_a2(pc, k, l)
```

## Arguments

- pc:

  The result of
  [`mvg_pieces`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
  with `derivs2`.

- k, l:

  Positions among the structure's free values.

## Value

A \\p \times p\\ numeric matrix.
