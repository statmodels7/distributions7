# Matrix Entries as the Default Interpretable Quantities

The distinct entries of the matrix
[`mv_sigma`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
returns, with a Jacobian from one central difference in each parameter.
This is what a family gets when it says nothing more specific: the
matrix on its own scale, named after the coordinates, rather than the
structure's coordinates.

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
