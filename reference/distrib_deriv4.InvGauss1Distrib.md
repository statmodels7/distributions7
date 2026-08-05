# Inverse-Gaussian Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Inverse-Gaussian log-density
(observed, or expected when `expected = TRUE`).

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
