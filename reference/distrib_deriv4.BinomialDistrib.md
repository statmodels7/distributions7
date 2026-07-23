# Binomial Analytical Fourth-Order Derivatives

Closed-form fourth-order derivative of the Binomial log-mass (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `BinomialDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivative.

## Value

A named list with the `mu_mu_mu_mu` component.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
