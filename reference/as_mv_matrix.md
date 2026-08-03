# Coerce a Multivariate Response to a Matrix

Puts `y` in the \\n \times p\\ form every multivariate method expects,
and refuses a response of the wrong width.

## Usage

``` r
as_mv_matrix(distrib, y)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- y:

  The response.

## Value

A numeric matrix with `distrib@n_dim` columns.

## Details

A plain vector of the right length is read as a single observation,
which is what makes `distrib_pdf(d, c(0, 0), theta)` mean what a reader
expects for a two-dimensional distribution. Anything else must already
be a matrix with one column per coordinate: guessing at the orientation
of an \\n \times p\\ matrix would silently transpose a square sample.
