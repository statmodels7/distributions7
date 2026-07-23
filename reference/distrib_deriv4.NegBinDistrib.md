# Negative Binomial Analytical Fourth-Order Derivatives

Closed-form fourth-order derivatives of the Negative Binomial log-mass.
The expected pure-\\\theta\\ derivative involves
\\\mathbb{E}\[\psi_3(Y+\theta)\]\\, evaluated by summation over the
support (with a far-tail correction).

## Arguments

- distrib:

  A `NegBinDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `theta`.

- expected:

  Logical; if `TRUE`, returns the expected fourth derivatives.

## Value

A named list of fourth-derivative component vectors.

## See also

[`negbin_distrib`](https://statmodels7.github.io/distributions7/reference/negbin_distrib.md)
