# Inverse-Gaussian Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Inverse-Gaussian log-density
(observed, or expected when `expected = TRUE`).

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `phi`.

- expected:

  Logical; if `TRUE`, returns the expected third derivatives.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list of third-derivative component vectors.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
