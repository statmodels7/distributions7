# Lognormal Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Lognormal log-density
(observed, or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `LognormalDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`lognormal_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
