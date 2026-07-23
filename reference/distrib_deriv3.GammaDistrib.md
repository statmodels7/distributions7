# Gamma Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Gamma log-density (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `GammaDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- expected:

  Logical; if `TRUE`, returns the expected third derivatives.

## Value

A named list of third-derivative component vectors.

## See also

[`gamma_distrib`](https://statmodels7.github.io/distributions7/reference/gamma_distrib.md)
