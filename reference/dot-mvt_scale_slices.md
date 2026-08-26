# Scale the Slices of an Array

Multiplies the \\i\\th slice of a \\p \times p \times n\\ array by
`v[i]`, which is how a term whose matrix part is the same at every
observation but whose scalar part is not enters an assembled derivative.
The multiplication is one vectorized recycle over the flattened array,
so nothing loops over the observations.

## Usage

``` r
.mvt_scale_slices(arr, v)
```

## Arguments

- arr:

  A \\p \times p \times n\\ numeric array.

- v:

  A numeric vector of length \\n\\. A shorter vector recycles silently,
  as the arithmetic underneath would.

## Value

A \\p \times p \times n\\ numeric array with the same dimensions as
`arr`.

## See also

[`.mvt_const_slices()`](https://statmodels7.github.io/distributions7/reference/dot-mvt_const_slices.md)
for the case where the matrix is what stays constant,
[`mv_outer_rows()`](https://statmodels7.github.io/distributions7/reference/mv_outer_rows.md)
for the array these are combined with, and
[`distrib_cross2_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md)
for the consumer.

## Examples

``` r
arr <- array(1, c(2, 2, 3))
distributions7:::.mvt_scale_slices(arr, c(1, 10, 100))[1, 1, ]
#> [1]   1  10 100

# Slice by slice it is the scalar times the slice.
a <- array(rnorm(12), c(2, 2, 3))
s <- distributions7:::.mvt_scale_slices(a, c(2, -1, 0.5))
all.equal(s[, , 2], -a[, , 2])
#> [1] TRUE
```
