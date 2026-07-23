# Laplace Probability Density Function

Computes the probability density function for the Laplace distribution:
\$\$f(y; \mu, b) = \dfrac{1}{2b}
\exp\left(-\dfrac{\|y-\mu\|}{b}\right)\$\$

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `b`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
