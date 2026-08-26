# Repeat a Constant Matrix Across Slices

Builds the \\p \times p \times n\\ array whose \\i\\th slice is
`v[i] * m`, the contribution of a term with a fixed matrix and a
per-observation scalar. It is the mirror of
[`.mvt_scale_slices()`](https://statmodels7.github.io/distributions7/reference/dot-mvt_scale_slices.md),
where the array varies and the scalar is applied to it; here the matrix
is the same everywhere and only the scalar moves.

## Usage

``` r
.mvt_const_slices(m, v)
```

## Arguments

- m:

  A \\p \times p\\ numeric matrix.

- v:

  A numeric vector of length \\n\\, which sets the third dimension of
  the result.

## Value

A \\p \times p \times n\\ numeric array, with `n = length(v)`.

## See also

[`.mvt_scale_slices()`](https://statmodels7.github.io/distributions7/reference/dot-mvt_scale_slices.md)
for the mirror case and
[`distrib_hess_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvStudentTDistrib.md)
for the consumer.

## Examples

``` r
m <- matrix(c(1, 2, 2, 4), 2, 2)
cs <- distributions7:::.mvt_const_slices(m, c(1, -1))
dim(cs)
#> [1] 2 2 2
cs[, , 2]
#>      [,1] [,2]
#> [1,]   -1   -2
#> [2,]   -2   -4

# Every slice is a multiple of the one matrix.
all.equal(cs[, , 1], m)
#> [1] TRUE
all.equal(cs[, , 2], -m)
#> [1] TRUE
```
