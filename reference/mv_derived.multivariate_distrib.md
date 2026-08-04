# Matrix Entries as the Default Interpretable Quantities

The distinct entries of the matrix
[`mv_sigma`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
returns, with a Jacobian from one central difference in each parameter.
This is what a family gets when it says nothing more specific: the
matrix on its own scale, named after the coordinates, rather than the
matrix parameter's coordinates.

The block a parameters7 family declares through
[`param_readable`](https://statmodels7.github.io/parameters7/reference/param_readable.html),
with its Jacobian widened from the free vector to the whole parameter
vector of the distribution.

## Usage

``` r
mv_param_block(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

A list as described in
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).

A list in the shape of
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
or `NULL`.

## Details

The free values of the structure occupy a contiguous stretch of
`distrib@params`, after the means and before anything the family adds of
its own, so widening the Jacobian is placing its columns in that stretch
and leaving the rest at zero: the quantities depend on no other
parameter. A family that declares nothing yields `NULL` and the summary
is what it was.
