# Negative Binomial Analytical Third-Order Derivatives

Closed-form third-order derivatives of the Negative Binomial log-mass.
The expected pure-\\\theta\\ derivative involves
\\\mathbb{E}\[\psi_2(Y+\theta)\]\\, evaluated by summation over the
support (with a far-tail correction).

## Arguments

- distrib:

  A `NegBin2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `theta`.

- expected:

  Logical; if `TRUE`, returns the expected third derivatives.

- threads:

  How many threads the kernel may use; below the measured internal
  threshold it stays sequential whatever the count says.

## Value

A named list of third-derivative component vectors.

## See also

[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
