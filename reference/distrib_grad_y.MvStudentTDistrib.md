# Multivariate Student t Response Gradient

\\\partial \ell / \partial y = -c\\\Sigma^{-1}(y-\mu)\\, the gaussian
expression with the family's weight in front of it.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

An \\n \times p\\ numeric matrix.
