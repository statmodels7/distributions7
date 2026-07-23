# Zero-Adjusted Continuous Analytical Observed Hessian

Observed Hessian of the zero-adjusted continuous model: the mixed blocks
are 0, the \\\pi\pi\\ block is \\-1/\pi^2\\ at \\y=0\\ and
\\-1/(1-\pi)^2\\ otherwise, and the \\\theta\theta\\ block is the
parent's Hessian at \\y \neq 0\\ (0 at \\y=0\\).

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list containing the vectors of second derivatives.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
