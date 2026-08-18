# Poisson Analytical Gradient

Computes the analytical gradient (first derivative) of the Poisson
log-probability with respect to the parameter \\\mu\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu}\$\$

## Arguments

- distrib:

  A `PoissonDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A list containing the vector of first derivatives.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
