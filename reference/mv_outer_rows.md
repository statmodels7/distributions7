# An Outer Product Per Observation

Builds the \\p \times p \times n\\ array whose \\i\\th slice is
\\a_ib_i^\top\\, the outer product of row \\i\\ of `a` with row \\i\\ of
`b`. The whole array is assembled by two column replications and one
permutation, with no loop over the observations.

## Usage

``` r
mv_outer_rows(a, b)
```

## Arguments

- a, b:

  Numeric matrices with the same dimensions, \\n\\ rows and \\p\\
  columns. The function does not validate them, and matrices of
  different dimensions signal the error of the arithmetic underneath.

## Value

A \\p \times p \times n\\ numeric array. It is symmetric slice by slice
only when `a` and `b` are the same matrix.

## See also

[`.mvt_scale_slices()`](https://statmodels7.github.io/distributions7/reference/dot-mvt_scale_slices.md)
for the other array helper and
[`distrib_cross2_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md)
for the consumer.

## Examples

``` r
a <- matrix(1:6, 3, 2)
b <- matrix(c(1, 0, -1, 2, 1, 0), 3, 2)
o <- distributions7:::mv_outer_rows(a, b)
dim(o)
#> [1] 2 2 3

# Slice i is the outer product of the two i-th rows.
o[, , 1]
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    4    8
outer(a[1, ], b[1, ])
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    4    8
all.equal(o[, , 2], outer(a[2, ], b[2, ]))
#> [1] TRUE

# If one matrix is given twice, every slice is symmetric.
s <- distributions7:::mv_outer_rows(a, a)
all(vapply(seq_len(dim(s)[3]),
           function(i) isSymmetric(s[, , i]), TRUE))
#> [1] TRUE
```
