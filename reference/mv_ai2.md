# The Second Derivative of the Inverse, by Position

As
[`mvg_a2()`](https://statmodels7.github.io/distributions7/reference/mvg_a2.md),
for the second derivative of \\\Sigma^{-1}\\ that
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md)
reads from the chart. Shared by the gaussian and the Student t.

## Usage

``` r
mv_ai2(pc, k, l)
```

## Arguments

- pc:

  Pieces carrying `ai2`, from
  [`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md)
  with `derivs2`.

- k, l:

  Positions among the structure's free values.

## Value

A \\p \times p\\ numeric matrix.
