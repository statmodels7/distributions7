# Transformed Third Derivatives

Exactly the parent's, evaluated at \\x = g^{-1}(y)\\: the Jacobian does
not depend on \\\theta\\, so it leaves every derivative in \\\theta\\
untouched.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list of the parent's parameters.

- expected:

  Logical; if `TRUE`, the expected derivatives.

## Value

A named list of derivative components.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
