# Gaussian Analytical Fourth-Order Derivatives

Computes the closed-form fourth-order partial derivatives of the
Gaussian log-density with respect to \\\mu\\ and \\\sigma\\ (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list of fourth-derivative component vectors.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
