# Default Mixed Derivatives for Continuous Distributions

Fallback: one central difference of
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
in each parameter (see
[`numerical_cross_y`](https://statmodels7.github.io/distributions7/reference/numerical_cross_y.md)).

## Arguments

- distrib:

  A `continuous_distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
