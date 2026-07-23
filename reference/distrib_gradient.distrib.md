# Default Numerical Gradient for `distrib` Objects

Fallback method: distributions that do not implement an analytical
gradient get one computed by central finite differences of the
log-density (see
[`numerical_gradient`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md)).

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list of gradient vectors.
