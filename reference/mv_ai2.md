# The Second Derivative of the Inverse, by Position

Reads the component of the second derivative array of \\\Sigma^{-1}\\
that
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md)
supplies, for the pair of free values \\(k, l)\\. The key is built from
the sorted pair of free names, so the result does not depend on the
order of the two positions. Shared by the gaussian and the Student t.

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

  Positions among the free values of the matrix parametrization.

## Value

A \\p \times p\\ numeric matrix.
