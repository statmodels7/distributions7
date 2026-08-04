# Transformed Probability Density Function

Change of variables: \\f_Y(y) = f_X(g^{-1}(y)) \cdot \|J(y)\|\\.
Computed in log space; singular log-densities are clamped to avoid
`Inf - Inf` during integration.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list of the parent's parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
