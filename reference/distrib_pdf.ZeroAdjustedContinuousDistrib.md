# Zero-Adjusted Continuous Probability Density Function

\$\$f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi) f_W(y;\theta) \\ \\ (y \neq
0)\$\$ (mixed density: point mass at 0 plus scaled continuous density).

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
