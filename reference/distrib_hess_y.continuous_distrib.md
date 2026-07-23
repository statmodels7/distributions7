# Default Response Hessian for Continuous Distributions

Fallback: \\\partial^2 \ell / \partial y^2\\ via finite differences (see
[`numerical_hess_y`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md)).

## Arguments

- distrib:

  A `continuous_distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A numeric vector.
