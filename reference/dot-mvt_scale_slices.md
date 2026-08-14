# Scale the Slices of an Array, and Repeat a Constant Matrix

`.mvt_scale_slices` multiplies the \\i\\th slice of an array by `v[i]`;
`.mvt_const_slices` returns the array whose \\i\\th slice is `v[i] * m`,
which is what a term constant in the observation contributes.

## Usage

``` r
.mvt_scale_slices(arr, v)

.mvt_const_slices(m, v)
```

## Arguments

- arr:

  A \\p \times p \times n\\ array.

- v:

  A numeric vector of length \\n\\.

- m:

  A \\p \times p\\ matrix.

## Value

A \\p \times p \times n\\ numeric array.
