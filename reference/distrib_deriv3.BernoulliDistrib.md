# Bernoulli Analytical Third-Order Derivatives

Closed-form third-order derivative of the Bernoulli log-mass (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- expected:

  Logical; if `TRUE`, returns the expected third derivative.

## Value

A named list with the `mu_mu_mu` component.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
