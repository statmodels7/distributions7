# Second-Response Mixed Derivatives of a Reparametrized Distribution

The parent's block carried by the first-order chain rule, exactly as
[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
is: a derivative in the response does not interact with a
reparametrization of the parameters.

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of the new parameters.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.
