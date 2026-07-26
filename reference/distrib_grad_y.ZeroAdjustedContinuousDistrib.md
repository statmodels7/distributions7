# Zero-Adjusted Continuous Response Gradient

\\\partial \ell / \partial y\\ equals the parent's for \\y \neq 0\\,
since the factor \\1-\pi\\ does not depend on \\y\\. At \\y = 0\\ the
log-density jumps to the atom and no derivative exists, so `NaN` is
returned.

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
