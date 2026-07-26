# Zero-Adjusted Continuous Response Hessian

\\\partial^2 \ell / \partial y^2\\ equals the parent's for \\y \neq 0\\
and is `NaN` at the atom.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A numeric vector.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
