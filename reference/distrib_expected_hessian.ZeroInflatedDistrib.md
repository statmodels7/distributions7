# Zero-Inflated Analytical Expected Hessian

Expected Hessian (negative Fisher information) of the zero-inflated
model, derived by decomposing the expectation over \\y=0\\ and \\y\>0\\
and using the parent's expected Hessian.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `zi`.

## Value

A list containing the vectors of expected second derivatives.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
