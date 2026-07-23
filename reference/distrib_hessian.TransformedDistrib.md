# Transformed Analytical Observed Hessian

The observed Hessian of the transformed model equals the parent's
observed Hessian evaluated at \\x = g^{-1}(y)\\.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list of the parent's parameters.

## Value

A list containing the vectors of second derivatives.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
