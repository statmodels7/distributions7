# The Pieces a Multivariate Gaussian Evaluates From

Assembles, once per call, the mean, the covariance, its inverse and its
log-determinant from a flat parameter vector, together with the matrix
parameter's derivative matrices when they are asked for.

## Usage

``` r
mvg_pieces(distrib, theta, derivs = FALSE, derivs2 = FALSE)
```

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

- derivs:

  Whether the matrix parameter's first derivative matrices are needed.

- derivs2:

  Whether its second derivatives are needed as well.

## Value

A list with `mu`, `sigma`, `sigma_inv`, `logdet`, `eta`, and optionally
`a` and `a2`, the derivatives of the covariance with respect to the free
values.

## Details

Whichever side the matrix parameter parametrises, the arithmetic below
is written in the covariance, so a precision structure is inverted once
here rather than at every use. The log-determinant follows the matrix
parameter's own, with its sign flipped for a precision, which is the one
place the two forms differ in anything but cost.
