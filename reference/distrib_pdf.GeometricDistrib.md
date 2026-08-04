# Geometric Probability Mass Function

\$\$P(Y = y; \mu) = \dfrac{1}{1+\mu}
\left(\dfrac{\mu}{1+\mu}\right)^{y}, \qquad y = 0, 1, 2, \dots\$\$

## Arguments

- distrib:

  A `GeometricDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
