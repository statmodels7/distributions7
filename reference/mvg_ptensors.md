# Precision Derivative Tensors of a Multivariate Gaussian

The precision's derivative tensors in the matrix parameter's free
values, orders 1 to 4, keyed by index tuple. For a precision
parametrization they are the parameter's own derivatives; for a
covariance they follow from repeated differentiation of the inverse, so
no expanded formula is transcribed and no term can be dropped. It takes
the pieces rather than the distribution, because the multivariate
Student t needs the same tensors of its scale matrix and there must be
one copy of the expansion: its first draft double counted the mixed
terms, which only a comparison against a stencil caught.

## Usage

``` r
mvg_ptensors(pc, order, inverted = FALSE)
```

## Arguments

- pc:

  The pieces, as returned by
  [`mvg_pieces`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
  or
  [`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md):
  a list carrying `s`, `eta` and `sigma_inv`.

- order:

  The highest order wanted.

- inverted:

  Whether the free values parametrize the precision, in which case the
  tensors are the parameter's own.

## Value

A list with the accessor `get`, the log-determinant sign and the pieces.
