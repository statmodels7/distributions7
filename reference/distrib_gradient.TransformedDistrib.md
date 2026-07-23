# Transformed Analytical Gradient

The Jacobian does not depend on the parameters, so the score of the
transformed model equals the parent's score evaluated at \\x =
g^{-1}(y)\\.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list of the parent's parameters.

## Value

A list containing the vectors of first derivatives.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
