# Zero-Adjusted Continuous Third Derivatives

There is no truncation constant, so away from the atom the \\\theta\\
derivatives are the parent's and the mixed ones vanish.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

## Value

A named list of derivative components.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
