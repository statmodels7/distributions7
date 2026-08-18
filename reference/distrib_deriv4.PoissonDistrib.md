# Poisson Analytical Fourth-Order Derivatives

Closed-form fourth-order derivative of the Poisson log-mass (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `PoissonDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivative.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list with the `mu_mu_mu_mu` component.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
