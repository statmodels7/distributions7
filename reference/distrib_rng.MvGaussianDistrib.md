# Multivariate Gaussian Generator

\\\mu + L z\\ with \\z\\ standard normal and \\LL^\top = \Sigma\\, the
factor taken from the matrix parameter where it parametrises the
covariance and from a factorisation of the inverse otherwise.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- n:

  The number of observations to draw.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

An \\n \times p\\ numeric matrix.
