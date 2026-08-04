# Zero-Adjusted Discrete Probability Mass Function

\$\$P(Y=0) = \pi, \qquad P(Y=y) =
(1-\pi)\dfrac{f(y;\theta)}{1-f(0;\theta)} \\ \\ (y\>0)\$\$

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `za`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of density values.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
