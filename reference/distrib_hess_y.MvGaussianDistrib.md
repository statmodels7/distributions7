# Multivariate Gaussian Response Hessian

\\\partial^2 \ell / \partial y \partial y^\top = -\Sigma^{-1}\\, the
same matrix at every observation, so it is returned once rather than
repeated.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

A \\p \times p\\ numeric matrix.
