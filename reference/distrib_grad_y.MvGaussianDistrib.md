# Multivariate Gaussian Response Gradient

\\\partial \ell / \partial y = -\Sigma^{-1}(y - \mu)\\, one row per
observation. The shape differs from the univariate case, where the
derivative in a scalar response is a vector: here it is an \\n \times
p\\ matrix.

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

An \\n \times p\\ numeric matrix.
