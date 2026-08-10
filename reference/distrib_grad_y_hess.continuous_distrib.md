# Default Second-Order Mixed Derivatives for Continuous Distributions

Fallback: one central difference of the analytic
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
or
[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
in each parameter (see
[`numerical_theta2_y`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)).

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

A named list keyed by parameter pair.
