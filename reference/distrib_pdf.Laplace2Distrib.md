# Laplace Density in Location and Rate

\$\$f(y; \mu, \lambda) = \dfrac{\lambda}{2}
\exp\left(-\lambda\|y-\mu\|\right)\$\$

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `lambda`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
