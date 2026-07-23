# Gamma Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Gamma log-density (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `GammaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
