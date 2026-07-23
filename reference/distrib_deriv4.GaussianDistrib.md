# Gaussian Analytical Fourth-Order Derivatives

Computes the closed-form fourth-order partial derivatives of the
Gaussian log-density with respect to \\\mu\\ and \\\sigma\\ (observed,
or expected when `expected = TRUE`).

## Arguments

- distrib:

  A `GaussianDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`gaussian_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md)
