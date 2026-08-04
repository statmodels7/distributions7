# Prefix a Structure's Free Names with the Matrix They Describe

Returns the structure's free names with `"sigma_"` or `"omega_"` in
front, according to which side of the model the structure parametrises.

## Usage

``` r
mv_prefixed_names(free_names, inverted = FALSE)
```

## Arguments

- free_names:

  The structure's free names.

- inverted:

  Whether the structure parametrises the precision.

## Value

A character vector.

## Details

The name of a free value says how the matrix is built, not which matrix
it is, so a covariance structure and a precision structure of the same
family produce identical names. They are different models — the inverse
of a compound-symmetry matrix is compound symmetry while the inverse of
an AR(1) is not AR(1) — and a printed table that does not distinguish
them leaves the reader to guess. The prefix is applied by the
distribution rather than by the structure, because the structure does
not know which side it has been handed to.
