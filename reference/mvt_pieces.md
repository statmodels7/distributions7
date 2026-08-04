# The Pieces a Multivariate t Evaluates From

Assembles the location, the scale matrix, its inverse, the
log-determinant and the degrees of freedom from a flat parameter vector,
with the structure's derivative matrices when they are needed.

## Usage

``` r
mvt_pieces(distrib, theta, derivs = FALSE, derivs2 = FALSE)
```

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

- derivs:

  Whether the first derivative matrices are needed.

- derivs2:

  Whether the second derivatives are needed as well.

## Value

A list with `mu`, `sigma`, `sigma_inv`, `logdet`, `nu`, `eta`, `p`, `s`,
and optionally `a` and `a2`.
