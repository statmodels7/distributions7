# Poisson Probability Mass Function

Computes the probability mass function for the Poisson distribution:
\$\$P(Y=y; \mu) = \dfrac{\mu^y e^{-\mu}}{y!}\$\$

## Arguments

- distrib:

  A `PoissonDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
