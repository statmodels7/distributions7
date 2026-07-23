# Transformed Analytical Expected Hessian

Since \\\ell_Y(\theta; y) = \ell_X(\theta; g^{-1}(y)) + \log\|J(y)\|\\
and the Jacobian does not depend on \\\theta\\, the expected Hessian of
the transformed model is *exactly* the parent's expected Hessian (the
expectation is just re-parameterized by the change of variables). No
Monte Carlo approximation is needed.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list of the parent's parameters.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
