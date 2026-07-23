# Default Response Gradient for Continuous Distributions

Fallback: \\\partial \ell / \partial y\\ via finite differences (see
[`numerical_grad_y`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md)).

## Arguments

- distrib:

  A `continuous_distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A numeric vector.
