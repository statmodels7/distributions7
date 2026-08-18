# Gamma Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Gamma log-density (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `Gamma2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list of fourth-derivative component vectors.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
