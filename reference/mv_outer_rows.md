# An Outer Product Per Observation

The \\p \times p \times n\\ array whose \\i\\th slice is
\\a_ib_i^\top\\, assembled without a loop over the observations.

## Usage

``` r
mv_outer_rows(a, b)
```

## Arguments

- a, b:

  Matrices of \\n\\ rows and \\p\\ columns.

## Value

A \\p \times p \times n\\ numeric array.
